# MCNetwork

This code uses a Kinetic Monte Carlo simulation to estimate the behaviour of variable range hopping in Dopant Networks. More information can be found [here](https://www.researchsquare.com/article/rs-757616/latest.pdf).

## Installation
This installation guide has been written using Ubuntu on WSL2 (Windows Subsystem for Linux), but should be usable on a range of linux operating systems.

Install the dependencies through apt (for Ubuntu/Debian):
```bash
sudo apt-get install cmake libboost-all-dev libhdf5-serial-dev
```

Run the following commands to build the project :
```bash
git clone https://github.com/MUTUEL/MCNetwork.git
git clone https://github.com/mfem/mfem.git
mkdir mfem-build && cd mfem-build
cmake ../mfem
make
mkdir ../MCNetwork/build && cd ../MCNetwork/build
env MFEM_DIR=../../mfem-build cmake ..
make
```

Set up your conda environment (assuming you have anaconda installed already):
```bash
conda create --name mcnetwork numpy scipy matplotlib h5py
conda activate mcnetwork
```

## Usage
Define the device parameters you want to use in `data/in.txt`. To make a new device with Monte Carlo optimization and random voltages go to the `data` folder in your shell and enter the following command:
```bash
../build/MCnetwork --mnd --optMC --rSV
```

Allowed options:

| option        | description                                                                                 |
|---------------|---------------------------------------------------------------------------------------------|
| --continue    | continues last optimization. optimization mode will be detected, no further options needed  |
| --mnd         | make new device                                                                             |
| --optMC       | optimize control voltages using Monte Carlo search                                          |
| --optGen      | optimize control voltages using genetic algorithm                                           |
| --optBasinHop | optimize control voltages using basin hopping                                               |
| --run         | just run control voltages defined in in.txt                                                 |
| --rSV         | random start voltages (only in combination with opt)                                        |
| --dir arg     | define working dir. has to contain 'in.txt'                                                 |
| --verbose     | additionally outputs occupation, swapps and time to 'additionalData.hdf5'                   |
| --help        | produce help message                                                                        |

The created device can be analyzed using the python scripts in the `scripts` folder.

## #TODO

### Technical:
- set "-D SWAPTRACKER" in cmake file
- improve argument parsing
- add log file
- support INT data format in datafile
- split up optimizer in child classes
- make one central hdf5 file to combine SWAPTRACKER and TiMETRACKER

### Optional
- count only output electrode current
- improve hash algorithm for storing mode
- in optimzer: store only control voltages
- in genetic algorithm: do not run unchanged genomes again
*/


## Contributors
- [Marlon Becker](https://github.com/MarlonBecker) [Westfälische Wilhelms-Universität Münster]
- [Jesse Bakker](https://github.com/Jesse-Bakker) [University of Twente]

## License
[APACHE LICENSE, VERSION 2.0](https://www.apache.org/licenses/LICENSE-2.0)