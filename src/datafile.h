#ifndef DATAFILE_H
#define DATAFILE_H
#include "H5Cpp.h"
#include <string>
#include <vector>

#include <iostream>
#include <boost/multi_array.hpp>
#include <iostream>
#include <functional>
#include <vector>
#include <map>
#include <stdio.h>
using namespace H5;

class DataFile
{

private:    
    std::map<std::string,int>    indexMap;
    std::string                  filename;    
    int                          numberOfDatasets=0;
    std::vector< hsize_t >       dims;
    std::vector< hsize_t * >     dimsf;
    std::vector< hsize_t * >     size;
    std::vector< hsize_t * >     offset;
    std::vector< hsize_t * >     chunk_dims;
    std::vector< hsize_t * >     maxdimsf;

    DataSpace spaceDummy;


public:
    DataFile(std::string _filename);
    void createDataset(std::string datasetName, std::vector<int> dimensions);
    void createAttribute(std::string, double);

    template<typename Data>
    void addData(std::string datasetName,Data data);

};


template<typename Data>
void DataFile::addData(std::string datasetName,Data data){
    int index=indexMap.at(datasetName);
    size  [index][0]++;
    offset[index][0]=size[index][0]-1;

    spaceDummy=DataSpace( dims[index], dimsf[index]);


    H5File * file= new H5File(filename, H5F_ACC_RDWR );

    DataSet dataset=file->openDataSet(datasetName);

    dataset.extend(size[index]);
    DataSpace fspace = dataset.getSpace ();
    fspace.selectHyperslab( H5S_SELECT_SET, dimsf[index], offset[index]);


    dataset.write( data , PredType::NATIVE_DOUBLE, spaceDummy, fspace );

    delete file;
}

#endif // DATAFILE_H
