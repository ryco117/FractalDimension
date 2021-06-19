﻿(*
FractalDimension - Experimental ray-marching based audio visualizer
Copyright (C) 2020  Ryan Andersen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
*)

module Program

open OpenTK.Windowing.Common // FOr the WindowBorder and WindowState oddly..
open OpenTK.Windowing.Desktop
open OpenTK.Graphics.OpenGL
open OpenTK.Mathematics
open OpenTK.Windowing.GraphicsLibraryFramework // For keyboard/mouse input

// Define type for storing a note
type Note = {freq: float32; mag: float32}
type NoteRange = Bass | Mids | High

// Fractal Types?!?!
type DistanceEstimate = None | Mandelbox | Mandelbulb | Klienian | Menger | IFS
let deToInt = function
| None -> 0
| Mandelbox -> 1
| Mandelbulb -> 2
| Klienian -> 3
| Menger -> 4
| IFS -> 5

let windowsSettings = NativeWindowSettings.Default
windowsSettings.NumberOfSamples <- 8
windowsSettings.Title <- "FractalDimension"

let toWorldSpace t noteType =
    let s = match noteType with
            | Bass -> System.Math.Pow(float t, 0.75)
            | Mids -> System.Math.Pow(float t, 0.55)
            | High -> System.Math.Pow(float t, 0.4)
    CubeFillingCurve.curveToCubeN 8 s

let icon =
    let imageData =
        let ico = new System.Drawing.Icon "FractalDimensionLogo.ico"
        use ms = new System.IO.MemoryStream ()
        let bitmap = ico.ToBitmap ()
        let bitmapData = bitmap.LockBits (System.Drawing.Rectangle (System.Drawing.Point(0, 0), bitmap.Size), System.Drawing.Imaging.ImageLockMode.ReadOnly, System.Drawing.Imaging.PixelFormat.Format32bppArgb)
        let scan = bitmapData.Scan0
        let stride = bitmapData.Stride
        for y = 0 to bitmapData.Height - 1 do
            for x = 0 to bitmapData.Width - 1 do
                let blue = System.Runtime.InteropServices.Marshal.ReadByte (scan, stride * y + 4 * x)
                let green = System.Runtime.InteropServices.Marshal.ReadByte (scan, stride * y + 4 * x + 1)
                let red = System.Runtime.InteropServices.Marshal.ReadByte (scan, stride * y + 4 * x + 2)
                let alpha = System.Runtime.InteropServices.Marshal.ReadByte (scan, stride * y + 4 * x + 3)
                ms.WriteByte red
                ms.WriteByte green
                ms.WriteByte blue
                ms.WriteByte alpha
        bitmap.UnlockBits bitmapData
        bitmapData.Width, bitmapData.Height, ms.ToArray ()
    Input.WindowIcon [|Input.Image imageData|]

type FractalDimension() =
    inherit GameWindow(GameWindowSettings.Default, windowsSettings, Icon = icon)

    // Window/Game World
    let mutable preFullscreenSize = Vector2i.One
    let mutable mouseHidden = false
    let mutable mouseLastMove = System.DateTime.UnixEpoch
    let cursorWait = 0.1
    // TODO: Remove this when OpenTK fixes MouseCursor.Empty
    let emptyCursor = new Input.MouseCursor(0, 0, 1, 1, Array.zeroCreate<byte> 4)
    let mutable tick = 0UL
    let mutable distanceEstimate = Mandelbox

    // Render
    let mutable renderShader = 0u
    let mutable vertexArrayObject = 0u

    // Demo
    let autoOrbitSpeed = 0.05f
    let orbitDist = 3.8f
    let mutable playTime = 0.
    let mutable kaleido = false
    let mutable kaleidoTime = -100.
    let mutable lastAngularChange = System.DateTime.UtcNow
    let mutable cubeAngularVelocity = Vector4(0.f, 1.f, 0.f, float32 autoOrbitSpeed)

    // Camera
    let mutable projectionInverse = Matrix4.Identity
    let mutable position = Vector3(0.f, 0.f, 2.21f)
    let mutable rotation = Quaternion(0.f, 0.f, 0.f, 1.f)
    let projectiveNear = 1.f
    let projectiveFar = 2.f

    // Audio
    let mutable audioResponsive = true
    let audioDisconnectCheckRate = 30UL
    let mutable complexZero = NAudio.Dsp.Complex ()
    do complexZero.X <- 0.f
    do complexZero.Y <- 0.f
    let noteZero = {freq = 0.f; mag = 0.f}
    let mutable previousBass = Array.create 4 [|complexZero|]
    let mutable previousBassIndex = 0
    let mutable smoothAverageBass = Vector3.Zero
    let mutable previousHigh = Array.create 8 noteZero
    let mutable previousHighIndex = 0
    let mutable smoothAverageHigh = Vector3.Zero
    let mutable smoothAverageHigh2 = Vector2.Zero
    let mutable previousMids = Array.create 8 noteZero
    let mutable previousMidsIndex = 0
    let mutable smoothAverageMids = Vector3.Zero
    let mutable smoothAverageMids2 = Vector2.Zero
    let mutable bassHoleTarget = Vector3(1.f, 0.f, 0.f)
    let mutable midsHoleTarget = Vector3(0.f, 1.f, 0.f)
    let mutable highHoleTarget = Vector3(0.f, 0.f, 1.f)
    let mutable lastVolume = 0.0001f
    let mutable bassHole = bassHoleTarget
    let mutable midsHole = midsHoleTarget
    let mutable highHole = highHoleTarget
    let onDataAvail samplingRate (complex: NAudio.Dsp.Complex[]) =
        if complex.Length > 0 then
            let mag (c: NAudio.Dsp.Complex) = sqrt(c.X*c.X + c.Y*c.Y)
            let freqResolution = samplingRate / float complex.Length
            let getStrongest maxCount delta (input: NAudio.Dsp.Complex[]) =
                let fLen = float32 input.Length
                let arr = Array.init input.Length (fun i -> {freq = (float32 i) / fLen; mag = mag input.[i]})
                let cmp {freq = _; mag = a} {freq = _; mag = b} = sign (b - a)
                let sorted = Array.sortWith cmp arr
                let rec getList acc size (arr: Note[]) =
                    if arr.Length = 0  || size = maxCount then
                        acc
                    else
                        let t = arr.[0].freq
                        let remaining, friends = Array.partition (fun {freq = s; mag = _} -> abs (t - s) > delta) arr
                        let m = Array.fold (fun acc {freq = _; mag = m} -> acc + m) 0.f friends
                        getList ({freq = t; mag = m}::acc) (size + 1) remaining
                List.toArray (List.rev (getList [] 0 sorted))
            let roundToInt f = int (round f)
            let bassStartFreq = 20.
            let bassEndFreq = 250.
            let midsStartFreq = 250.
            let midsEndFreq = 3000.
            let highStartFreq = 3000.
            let highEndFreq = 15000.
            let bassStart = roundToInt (bassStartFreq / freqResolution)
            let bassEnd = roundToInt (bassEndFreq / freqResolution)
            let midsStart = roundToInt (midsStartFreq / freqResolution)
            let midsEnd = roundToInt (midsEndFreq / freqResolution)
            let highStart = roundToInt (highStartFreq / freqResolution)
            let highEnd = roundToInt (highEndFreq / freqResolution)
            let bassArray = Array.sub complex bassStart (bassEnd - bassStart)
            let bassNotes = getStrongest 1 0.4f bassArray
            let midsNotes = getStrongest 1 0.2f (Array.sub complex midsStart (midsEnd - midsStart))
            let highNotes = getStrongest 1 0.3f (Array.sub complex highStart (highEnd - highStart))
            let avgLastBassMag x =
                let mutable s = 0.f
                for i = 0 to previousBass.Length - 1 do
                    s <- s +
                        if previousBass.[i].Length = 0 then
                            0.f
                        else
                            let j =
                                let j = int (round (x * float32 previousBass.[i].Length))
                                if j >= previousBass.[i].Length then previousBass.[i].Length - 1 else j
                            mag previousBass.[i].[j]
                s / float32 previousBass.Length
            let volume = 
                let summer a = Array.sumBy (fun n -> n.mag) a
                summer bassNotes + summer midsNotes + summer highNotes
            let mutable canJerk = true
            let minimumBassForJerk = 0.05f
            let autoOrbitJerk = 0.185f
            let minimumBass = 0.0075f
            for i = 0 to bassNotes.Length - 1 do
                if bassNotes.[i].mag > minimumBass && bassNotes.[i].mag > 1.25f * avgLastBassMag bassNotes.[i].freq then
                    bassHoleTarget <- (1.f - minimumBass / bassNotes.[i].mag) * (toWorldSpace bassNotes.[i].freq Bass)
                if canJerk &&
                    bassNotes.[i].mag > minimumBassForJerk &&
                    (let t = (System.DateTime.UtcNow - lastAngularChange).TotalSeconds in t > 0.66 / float bassNotes.[i].mag || t > 1.5) &&
                    bassNotes.[i].mag > 10.f * avgLastBassMag bassNotes.[i].freq then
                    cubeAngularVelocity <- Vector4(
                        (toWorldSpace bassNotes.[i].freq Bass).Normalized(),
                        (sqrt volume) * autoOrbitJerk)
                    lastAngularChange <- System.DateTime.UtcNow
                    canJerk <- false
            let minimumMids = 0.001f
            for i = 0 to midsNotes.Length - 1 do
                if midsNotes.[i].mag > minimumMids then
                    let worldSpace = toWorldSpace midsNotes.[i].freq Mids
                    midsHoleTarget <- (1.f - minimumMids / midsNotes.[i].mag) * worldSpace
                    previousMids.[previousMidsIndex] <- {freq = System.MathF.Pow(midsNotes.[i].freq, 0.5f); mag = midsNotes.[i].mag}
                    previousMidsIndex <- (previousMidsIndex + 1) % previousMids.Length
            let minimumHigh = 0.0005f
            for i = 0 to highNotes.Length - 1 do
                if highNotes.[i].mag > minimumHigh then
                    let worldSpace = toWorldSpace highNotes.[i].freq High
                    highHoleTarget <- (1.f - minimumHigh / highNotes.[i].mag) * worldSpace
                    previousHigh.[previousHighIndex] <- {freq = System.MathF.Pow(highNotes.[i].freq, 0.75f); mag = 1.775f * highNotes.[i].mag}
                    previousHighIndex <- (previousHighIndex + 1) % previousHigh.Length
            lastVolume <- volume
            previousBass.[previousBassIndex] <- bassArray
            previousBassIndex <- (previousBassIndex + 1) % previousBass.Length
    let audioOutCapture = new EzSound.AudioOutStreamer(onDataAvail, fun () -> lastVolume <- 0.0001f)
    override this.OnLoad () =
        // Set default values
        GL.Hint(HintTarget.MultisampleFilterHintNv, HintMode.Nicest)
        GL.ClearColor (0.f, 0.f, 0.f, 1.f)
        GL.PolygonMode(MaterialFace.Front, PolygonMode.Fill)
        GL.DepthFunc DepthFunction.Lequal
        GL.Enable EnableCap.DepthTest
        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha)
        GL.Enable EnableCap.Blend

        // Create fractal window
        renderShader <- EzShader.CreateShaderProgram "vert.glsl" "frag.glsl"
        vertexArrayObject <- GL.GenVertexArray ()
        GL.BindVertexArray vertexArrayObject

        let vertexBufferObject = GL.GenBuffer ()
        GL.BindBuffer(BufferTargetARB.ArrayBuffer, vertexBufferObject)
        let vertices = [|
            -1.0f;  1.0f;
             1.0f;  1.0f;
            -1.0f; -1.0f;
             1.0f; -1.0f
        |]
        GL.BufferData(BufferTargetARB.ArrayBuffer, vertices, BufferUsageARB.StaticDraw)

        GL.UseProgram renderShader
        GL.VertexAttribPointer(0u, 2, VertexAttribPointerType.Float, false, 0, nativeint 0)
        GL.EnableVertexAttribArray 0u
        GL.Uniform1f(GL.GetUniformLocation(renderShader, "mandelboxScale"), -1.5f)
    override this.OnRenderFrame eventArgs =
        let deltaTime = eventArgs.Time
        GL.Clear (ClearBufferMask.ColorBufferBit ||| ClearBufferMask.DepthBufferBit)

        // Update the rotatiom amd angular momentum of the camera
        do
            playTime <- playTime + System.Math.Pow(float lastVolume, 0.65) * deltaTime
            let deltaTime = float32 deltaTime
            let w = cubeAngularVelocity.W
            let theta = w * deltaTime
            let r = Quaternion(sin theta * cubeAngularVelocity.Xyz, cos theta)
            rotation <- (rotation * r).Normalized()
            cubeAngularVelocity <- Vector4(cubeAngularVelocity.Xyz, w + (autoOrbitSpeed - w) * (1.f - exp (-deltaTime/2.75f)))

        // Update view matrix 
        let mutable viewMatrix = projectionInverse * Matrix4.CreateFromQuaternion rotation

        // Update position of note-vectors
        let smooth = (1.f - exp (float32 deltaTime / -1.4f))
        bassHole <- bassHole + (bassHoleTarget - bassHole) * smooth
        midsHole <- midsHole + (midsHoleTarget - midsHole) * smooth
        highHole <- highHole + (highHoleTarget - highHole) * smooth

        // Update the rotation amd angular momentum of the camera
        do
            let modulo = 15.f * MathHelper.Pi
            let half = modulo / 2.f
            let m = abs (((float32 playTime + half) % modulo) / half - 1.f)
            let d = orbitDist - 2.1f * m
            let projectionConstant = d*(projectiveFar+projectiveNear)/(projectiveFar-projectiveNear) - (2.f*projectiveFar*projectiveNear)/(projectiveFar-projectiveNear)
            position <- -(Vector4(0.f, 0.f, projectionConstant, d) * viewMatrix).Xyz

        // Perfrom render
        GL.BindVertexArray vertexArrayObject
        GL.UseProgram renderShader

        // Set the shader uniforms
        let viewMatrixBytes = Array.init 16 (fun i -> viewMatrix.Item (i / 4, i % 4))
        GL.UniformMatrix4f(GL.GetUniformLocation(renderShader, "projectiveInverse"), 1, false, viewMatrixBytes)
        GL.Uniform3f(GL.GetUniformLocation(renderShader, "cameraPosition"), &position)
        GL.Uniform3f(GL.GetUniformLocation(renderShader, "bassHole"), &bassHole)
        GL.Uniform3f(GL.GetUniformLocation(renderShader, "midsHole"), &midsHole)
        GL.Uniform3f(GL.GetUniformLocation(renderShader, "highHole"), &highHole)

        GL.Uniform1i(GL.GetUniformLocation(renderShader, "deType"), deToInt distanceEstimate)
        if distanceEstimate <> None then
            let magic = if distanceEstimate = Mandelbulb then -1. else -cos (playTime / 10.)
            GL.Uniform1f(GL.GetUniformLocation(renderShader, "magicNumber"), float32 magic)
            let magicLin = (playTime / 20.) % 2.
            GL.Uniform1f(GL.GetUniformLocation(renderShader, "magicNumberLin"), float32 (min magicLin (2. - magicLin)))
            let scale = -(0.15 * (asin (-cos (playTime / 6.)) + 1.) + 1.95)
            GL.Uniform1f(GL.GetUniformLocation(renderShader, "mandelboxScale"), float32 scale)
            let kaleidoscope =
                let t = min ((playTime - kaleidoTime) / 0.6) 1.
                sqrt (if kaleido then t else (1. - t))
            GL.Uniform1f(GL.GetUniformLocation(renderShader, "kaleido"), float32 kaleidoscope)

        let smoothScale (v: Vector3) (arr: Note[]) r =
            let smooth = (1.f - exp (float32 deltaTime / -6.f))
            let mutable avg = Vector3.Zero
            for v in arr do
                avg <- avg + toWorldSpace v.freq r
            avg <- avg / (float32 arr.Length)
            v + (avg - v) * smooth
        smoothAverageHigh <- smoothScale smoothAverageHigh previousHigh NoteRange.High
        GL.Uniform3f(GL.GetUniformLocation(renderShader, "avgHighHole"), &smoothAverageHigh)
        smoothAverageMids <- smoothScale smoothAverageMids previousMids NoteRange.Mids
        GL.Uniform3f(GL.GetUniformLocation(renderShader, "avgMidsHole"), &smoothAverageMids)

        let smoothScale2 (v: Vector2) (arr: Note[]) =
            let mutable avg = Vector2.Zero
            for v in arr do
                let theta = MathHelper.TwoPi * v.freq
                avg <- avg + System.MathF.Pow(v.mag, 0.425f) * Vector2(System.MathF.Cos theta, System.MathF.Sin theta)
            avg <- avg / (float32 arr.Length)
            v + (avg - v) * smooth
        smoothAverageHigh2 <- smoothScale2 smoothAverageHigh2 previousHigh
        GL.Uniform2f(GL.GetUniformLocation(renderShader, "avgHighPlot"), &smoothAverageHigh2)
        smoothAverageMids2 <- smoothScale2 smoothAverageMids2 previousMids
        GL.Uniform2f(GL.GetUniformLocation(renderShader, "avgMidsPlot"), &smoothAverageMids2)

        GL.DrawArrays(PrimitiveType.TriangleStrip, 0, 4)
        this.Context.SwapBuffers ()

        // Occasionally check if audio out stopped
        if audioResponsive && tick % audioDisconnectCheckRate = 0UL && audioOutCapture.Stopped () then
            audioOutCapture.StartCapturing ()

        if (System.DateTime.UtcNow.Subtract mouseLastMove).TotalSeconds > cursorWait then
            if not mouseHidden then
                this.Cursor <- emptyCursor
                mouseHidden <- true
        else
            if mouseHidden then
                this.Cursor <- Input.MouseCursor.Default
                mouseHidden <- false

        base.OnRenderFrame eventArgs
    override this.OnKeyDown e =
        match e.Key, (e.Alt, e.Shift, e.Control), e.IsRepeat with
        // Quit
        | Keys.F4, (true, false, false), _
        | Keys.Escape, _, _ -> exit 0
        // Fullscreen
        | Keys.Enter, (true, false, false), false
        | Keys.F11, (false, false, false), false ->
            if this.WindowBorder = WindowBorder.Hidden then
                this.WindowState <- WindowState.Normal
                this.WindowBorder <- WindowBorder.Resizable
                this.Size <- preFullscreenSize
            else
                preFullscreenSize <- this.Size
                this.WindowBorder <- WindowBorder.Hidden
                this.WindowState <- WindowState.Fullscreen
        | Keys.R, (false, false, false), false ->
            if audioResponsive then
                audioResponsive <- false
                if audioOutCapture.Capturing () then
                    audioOutCapture.StopCapturing ()
            else
                audioResponsive <- true
                audioOutCapture.Reset ()
        | Keys.D0, _, false ->
            distanceEstimate <- None
        | Keys.D1, _, false ->
                distanceEstimate <- Mandelbox
        | Keys.D2, _, false ->
            distanceEstimate <- Mandelbulb
        | Keys.D3, _, false ->
            distanceEstimate <- Klienian
        | Keys.D4, _, false ->
            distanceEstimate <- Menger
        | Keys.D5, _, false ->
            distanceEstimate <- IFS
        | Keys.Space, _, false ->
            kaleido <- not kaleido
            kaleidoTime <- playTime
        | _ -> base.OnKeyDown e
    override _.OnMouseMove move =
        mouseLastMove <- System.DateTime.UtcNow
        base.OnMouseMove move
    override _.OnResize size =
        GL.Viewport (0, 0, size.Width, size.Height)
        projectionInverse <- Matrix4.CreatePerspectiveFieldOfView(MathHelper.PiOver2, float32 size.Width / float32 size.Height, projectiveNear, projectiveFar).Inverted()
        base.OnResize size
    override _.OnClosing _ = exit 0

[<EntryPoint>]
let main _ =
    use world = new FractalDimension()
    world.Run ()
    0