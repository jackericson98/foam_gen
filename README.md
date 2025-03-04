# foam_gen

#### John Ericson
#### Georgia State University
#### 2025

#

Welcome to `foam_gen`, an interactive simulated foam generator. This project is designed to create 3-dimensional random ensembles of spheres, mimicking the bubbles in foam based on user-defined parameters. This program was initially designed to be used 

## Features

- Generate 3D random ensembles of spheres.
- Customize the size, distribution, and interaction of spheres.
- Export the generated foam structure for further analysis or visualization.

## Installation

To get started with `foam_gen`, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/jackericson98/foam_gen.git
   ```
2. Navigate to the project directory:
   ```bash
   cd foam_gen
   ```
4. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```
   or
   
   ```bash
   python -m pip install -r requirements.txt
   ```

## Usage

To run `foam_gen` type the following into the terminal from the `foam_gen` directory:

```bash
python foam_gen.py
```
This will trigger the foam_gen GUI:

![GUI_Image](https://github.com/user-attachments/assets/b31fb1ec-0f2f-4fd3-8af5-a1d4a0a170b6)

The GUI can then be used to define your foam and outputs. When ready select the `Create Foam` button and the foam will be created under the set paramaters. These paramaters are described below. 

## Parameters

1. Average: The mean radius of the set of balls. 
2. Coefficient of Variation: The poly dispersity of the balls. This has different meanings based on the different distributions, but generally refers to how wide the distribution is. 
3. Number of Balls: How many output balls are in the final set. (Note if the periodic boundary condition is true, the numbe of balls in the output file will be 27x this)
4. Density: How packed the balls are in the output file. This refers to the amount of free space in the retaining box vs how much volume the balls occupy. 
5. Overlap: How much the balls are allowed to cross each other. This parameter is in units of smaller ball radius. If two balls are placed near eachother the smaller radius ball is used as a guide for how close they can be. For example, if ball 1 has a radius of 2.5 and ball 2 has a radius of 0.5 and the overlap value is set to 1, the balls centers must be at least 2.0 away from each other. 
6. Distribution: This determines how the radii are distributed. Given the mean ball size and the CV the chosen distribution can be manipulated. 
7. Periodic Boundary: Whether the balls can interact with the balls on the opposite wall (think pacman). This affects the placement of the balls. To be able to better visualize how the balls interact with their opposing wall balls, the program outputs balls that surround the main set with the same orientatation and radii, just moved by one box length over. 
8. Atandard Atomic Radii: If the user wants only radii that correspond to atomic Van der Waal's radii (e.g. C = 1.6, H = 1.2), this option should be selected. foam_gen will find the closest radii from the set of radii in the radii file to a set ball's radius in the distribution specified and then give that ball the corresponding element.

## Outputs

There are 4 major outputs from a `foam_gen` run:

1. A `.txt` file with the balls location and radii. This also holds information about the specific atoms in the grouping if the 

## Contributing


## License

MIT License

Copyright (c) 2022 John Ericson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

   
