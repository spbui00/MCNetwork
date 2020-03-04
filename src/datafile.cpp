#include "datafile.h"


DataFile::DataFile(std::string _filename){

    filename = _filename;
    
    H5File * file= new H5File( filename, H5F_ACC_TRUNC );

    delete file;
}

void DataFile::createDataset(std::string datasetName, std::vector<int> dimensions){

    H5File * file= new H5File(filename, H5F_ACC_RDWR );

    indexMap[datasetName]=numberOfDatasets;
    numberOfDatasets++;
    dims      .push_back(dimensions.size()+1);
    dimsf     .push_back(new hsize_t[dims[numberOfDatasets-1]]);
    size      .push_back(new hsize_t[dims[numberOfDatasets-1]]);
    offset    .push_back(new hsize_t[dims[numberOfDatasets-1]]);
    chunk_dims.push_back(new hsize_t[dims[numberOfDatasets-1]]);
    maxdimsf  .push_back(new hsize_t[dims[numberOfDatasets-1]]);
    
    dimsf     [numberOfDatasets-1][0]=1;
    size      [numberOfDatasets-1][0]=0;
    offset    [numberOfDatasets-1][0]=0;
    chunk_dims[numberOfDatasets-1][0]=1;
    maxdimsf  [numberOfDatasets-1][0]=H5S_UNLIMITED;

    for (size_t i = 1; i < dims[numberOfDatasets-1]; i++){
        dimsf     [numberOfDatasets-1][i]=dimensions[i-1];
        size      [numberOfDatasets-1][i]=dimensions[i-1];
        offset    [numberOfDatasets-1][i]=0;
        chunk_dims[numberOfDatasets-1][i]=dimensions[i-1];
        maxdimsf  [numberOfDatasets-1][i]=H5S_UNLIMITED;
    }

    DataSpace dataspace( dims[numberOfDatasets-1], size[numberOfDatasets-1], maxdimsf[numberOfDatasets-1]);  


    DSetCreatPropList cparms;
    cparms.setChunk( dims[numberOfDatasets-1], chunk_dims[numberOfDatasets-1] );
    int fill_val = 0;
    cparms.setFillValue( PredType::NATIVE_DOUBLE, &fill_val);
   
    DataSet dataset = file->createDataSet( datasetName, PredType::NATIVE_DOUBLE, dataspace,cparms);
    
    delete file;
}



void DataFile::createAttribute(std::string attrName, double val)
{
    H5File file= H5File(filename, H5F_ACC_RDWR );

    const hsize_t dims=1;
    DataSpace* dspace = new DataSpace(1,&dims);
    Attribute attr = file.createAttribute(attrName, PredType::NATIVE_DOUBLE, *dspace);
    delete dspace;
    
    attr.write(PredType::NATIVE_DOUBLE,& val);
}

