![FractalDimension](TheFractalDimension/FractalDimensionLogo.ico)
# FractalDimension

### About the project
This project is derived from my other experimental audio visualizer, [ColouredSugar](https://github.com/ryco117/ColouredSugar), and owes much of its source code to that project.

An experimental audio visualizer written in F# and using the library [OpenTK](https://github.com/opentk/opentk) as a wrapper for the OpenGL API.
The fractal is rendered using a technique called [ray-marching](http://blog.hvidtfeldts.net/index.php/2011/06/distance-estimated-3d-fractals-part-i/).
I chose the open source library [NAudio](https://github.com/naudio/NAudio) for handling the audio streams, despite the fact that OpenTK comes with a wrapper for OpenAL, 
because I could not determine a way to stream the current audio out as a capture in OpenAL.
NAudio also provides the Fast Fourier Transform I used.

### Visualizer Technical Description
FractalDimension uses the [WASAPI loopback](https://docs.microsoft.com/en-us/windows/win32/coreaudio/loopback-recording) 
feature to capture the current audio out stream. The audio is then read in chunks as they become available 
and a fast fourier transform is applied to the sampled audio stream. The resulting frequency domain is partitioned into three unequal parts, 
`bass`, `mids`, and `highs`. For each of these sub domains, the discrete frequency bins with the largest magnitudes are determined 
(ie. the tones most present in the sound wave are determined). 
Then, based on the frequency and magnitude of the detected notes, various parameters of the environment are altered.

### Precompiled Executable (Windows)
- Go to the [releases page](https://github.com/ryco117/FractalDimension/releases) and determine the latest release/version
- Download the `FractalDimension-precompiled-Win64.zip` archive file
- Extract the folder within
- Run the `FractalDimension.exe` executable within

### Example Videos (Click for links)
##### Bassnectar Demo
[![Bassnectar Demo](https://img.youtube.com/vi/OQv0IDQd8H8/0.jpg)](https://www.youtube.com/watch?v=OQv0IDQd8H8 "Bassnectar Demo")
##### Queen Demo
[![Queen Demo](https://img.youtube.com/vi/TabDKU24C6s/0.jpg)](https://www.youtube.com/watch?v=TabDKU24C6s "Queen Demo")
##### Styx Demo
[![Styx Demo](https://img.youtube.com/vi/Eo2TI3FJtoM/0.jpg)](https://www.youtube.com/watch?v=Eo2TI3FJtoM "Styx Demo")
##### Luigi Sambuy Demo
[![L.Dre Demo](https://img.youtube.com/vi/_9rtE6Lh7aY/0.jpg)](https://www.youtube.com/watch?v=_9rtE6Lh7aY "Luigi Sambuy Demo")