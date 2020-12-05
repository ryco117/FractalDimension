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
##### Death Cab for Cutie
[![Death Cab for Cutie Demo](https://img.youtube.com/vi/p_YMFBeia4w/0.jpg)](https://www.youtube.com/watch?v=p_YMFBeia4w "Death Cab for Cutie Demo")
##### Madeon Demo
[![Madeon Demo](https://img.youtube.com/vi/f5pirPzv0Yk/0.jpg)](https://www.youtube.com/watch?v=f5pirPzv0Yk&t=2 "Madeon Demo")
##### Santana Demo
[![Santana Demo](https://img.youtube.com/vi/yuSlCbhd97U/0.jpg)](https://www.youtube.com/watch?v=yuSlCbhd97U&t=4 "Santana Demo")
##### L.Dre Demo
[![L.Dre Demo](https://img.youtube.com/vi/FVB3HRSrhi8/0.jpg)](https://www.youtube.com/watch?v=FVB3HRSrhi8&t=2 "L.Dre Demo")