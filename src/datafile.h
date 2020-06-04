#ifndef DATAFILE_H
#define DATAFILE_H
#include "H5Cpp.h"
#include <chrono>
#include <string>
#include <thread>
#include <vector>

#include <boost/multi_array.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <stdio.h>
#include <vector>
using namespace H5;

/*!
    container class for hdf5 output file. 
 */
class DataFile {
public:
    DataFile(std::string filename, bool createNew);
    void createDataset(std::string datasetName, std::vector<int> dimensions);
    void createAttribute(std::string, double);

    template <typename Data>
    void addData(std::string datasetName, Data data);
    void shrinkDataset(std::string datasetName, int dataPointsToDelete);

    std::unique_ptr<std::vector<double>> readFullDataset(std::string datasetName);
    std::unique_ptr<std::vector<double>> readDatasetSlice(std::string datasetName, int index);

    bool checkDataSetExists(std::string datasetName);

private:
    static std::map<std::string, int> indexMap; /*!< maps datset names to index to find properties, e.g. dims \n has to be static to be accessible by  C API in lambda function in constructor DataFile() */
    std::string filename;
    int numberOfDatasets = 0;
    std::vector<hsize_t> dims; /*!< see hdf5 c++ examples */
    std::vector<hsize_t*> dimsf; /*!< see hdf5 c++ examples */
    std::vector<hsize_t*> size; /*!< see hdf5 c++ examples */
    std::vector<hsize_t*> offset; /*!< see hdf5 c++ examples */
    std::vector<hsize_t*> chunk_dims; /*!< see hdf5 c++ examples */
    std::vector<hsize_t*> maxdimsf; /*!< see hdf5 c++ examples */
};

/*!
    extends existing dataset by one slice in dimension 0.
    \param  datasetName dataset must be created beforehand 
    \param data c-array !!! currently ONLY DOUBLE array supported !! i.e. ints will raise no error, but produce incorrect data!
 */
template <typename Data>
void DataFile::addData(std::string datasetName, Data data)
{
    int index = indexMap.at(datasetName);
    size[index][0]++;
    offset[index][0] = size[index][0] - 1;

    H5File file;

    int waitTime = 0;
tryAgian:;
    try {
        file = H5File(filename, H5F_ACC_RDWR);
    } catch (FileIException const& e) {
        waitTime++;
        if (waitTime == 1) {
            std::cout << "------------> Unable to open file <------------" << std::endl;
            std::cout << "---------------> error message: <--------------" << std::endl;
            e.clearErrorStack();
            std::cout << "--------> trying again in 1 second <-----------" << std::endl;
        } else {
            std::cout << "--------> already waiting for " << waitTime << " seconds <----------" << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::seconds(1));
        goto tryAgian;
    }

    //get data space
    DataSet dataset = file.openDataSet(datasetName);
    dataset.extend(size[index]);
    DataSpace fspace = dataset.getSpace();
    fspace.selectHyperslab(H5S_SELECT_SET, dimsf[index], offset[index]);

    //get mem space
    DataSpace memSpace(dims[index], dimsf[index]);

    dataset.write(data, PredType::NATIVE_DOUBLE, memSpace, fspace);
}

#endif // DATAFILE_H
