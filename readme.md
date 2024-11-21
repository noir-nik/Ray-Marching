## Distance-aided Ray Marching
### Rendering of implicit surfaces with GPU acceleration

## Example scene
![](/images/out.jpg)

## Features
  * Anti-aliasing (super sampling x4)
  * Sharp shadows
  * Reflections
  * Ambient Occlusion
  * Fractal surfaces:
	* Menger sponge
    * 3D Sierpinski triangle
    * Mandelbulb
  * Use of Constructive Solid Geometry
  * Unique color for each primitive
  * Smooth operations

## Building

```sh
make
```

## Running

```sh
make run
```

Rendered images in BMP format should appear in the `output` folder

## Render times
Measurements are done with i5 10400f (OpenMP) and RTX 3080 (OpenGL):

### No anti-aliasing:
```sh
CPU render - 1 thread:   42576.211 ms
CPU render - 12 threads:  5107.263 ms
GPU Render time:            13.889 ms
Copy + GPU render time:     14.385 ms
Copy time:                   0.496 ms
```
### AA x4:
```sh
CPU render - 1 thread:   168946.344 ms
CPU render - 12 threads:  20155.684 ms
GPU Render time:             48.111 ms
Copy + GPU render time:      48.608 ms
Copy time:                    0.497 ms
```