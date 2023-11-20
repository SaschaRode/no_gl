# noGL - software renderer

This software renderer was developed in two weeks and the goal of the project was to find a design that could be easily ported to GPU.

There is also a small example viewer that comes with the project.

Currently only Windows and MSVC is supported.

## What are the Features?

- Depth testing
- Alpha blending
- Perspective-correct interpolation
- Backface culling
- Clipping
- Fragment shader support
- Blinn-Phong lighting
- Normal mapping

## How to Build locally?

Start a Visual Studio x64 Command Prompt and navigate to the project's root directory.
Execute the build.bat file and you should be ready to go.

## Credits

- stb_image: Sean Barret (nothings) - https://github.com/nothings/stb
- objzero: Jonathan Young - https://github.com/jpcy/objzero