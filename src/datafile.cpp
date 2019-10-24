#include "datafile.h"



using namespace H5;


// DataFile::~DataFile()
// {
//     flush();
// }


// DataFile::DataFile(Lipidsystem& lipidsystem,std::shared_ptr<InputFile> _inputfile)
// {
//     inputfile=_inputfile;
    
//     getterMap["orderPara"]=std::bind(&Lipidsystem::getOrderParas,&lipidsystem);
//     getterMap["Type"]=std::bind(&Lipidsystem::getTypes,&lipidsystem);
//     getterMap["IDs"]=std::bind(&Lipidsystem::getIDs,&lipidsystem);


//     NX=inputfile->width;
//     NY=inputfile->height;
//     maxBufferSize=1024*1024*inputfile->paras.at("maxBufferSize");

// }


// void DataFile::createFile()
// {
//    H5File file( filename, H5F_ACC_TRUNC );
    
//     dimsf[0] = 1;
//     dimsf[1] = NX;
//     dimsf[2] = NY;

//     size[0] = 0;
//     size[1] = NX;
//     size[2] = NY;

//     chunk_dims[0] = 1;
//     chunk_dims[1] = NX;
//     chunk_dims[2] = NY;


//     for(auto& o: inputfile->outs)
//     {
//         createDataset(o,file);
//     }
    
//     for (auto const& mapitem: inputfile->paras) createAttribute(mapitem.first,mapitem.second,file);
    

// }

// void DataFile::createDataset(std::string datasetName, H5File& file)
// {
//     DataSpace dataspace( 3, size, maxdimsf);  

      
//     DSetCreatPropList cparms;
//     cparms.setChunk( 3, chunk_dims );
//     int fill_val = 0;
//     cparms.setFillValue( PredType::NATIVE_INT, &fill_val);
   
    
//     DataSet dataset = file.createDataSet( datasetName, PredType::NATIVE_INT, dataspace,cparms);
    
//     buffer.push_back(boost::multi_array<int,3>(boost::extents[0][NX][NY]));
// }

// void DataFile::createAttribute(std::string attrName, double val, H5File& file)
// {
//     double* data;
//     data=&val;
//     const hsize_t dims=1;
//     DataSpace* dspace = new DataSpace(1,&dims);
//     Attribute attr = file.createAttribute(attrName, PredType::NATIVE_DOUBLE, *dspace);
//     delete dspace;
    
//     attr.write(PredType::NATIVE_DOUBLE,data);
    
// }

// void DataFile::writeStep()
// {
//     #ifndef NDEBUG
//     std::cout<<"DataFile::writeStep"<<std::endl;
//     #endif

//     bufferLen++;
//     for(int i=0; i<inputfile->outs.size();i++)
//     {
//         buffer[i].resize(boost::extents[bufferLen][NX][NY]);
//         buffer[i][bufferLen-1]=getterMap.at(inputfile->outs[i])();

//     }
//     bufferSize+=NX*NY*sizeof(int)*inputfile->outs.size();
//     if(bufferSize>=maxBufferSize) flush();

// }

// void DataFile::flush()
// {

//     H5File file( filename, H5F_ACC_RDWR );
    
//     images+=bufferLen;
//     dimsf[0]=bufferLen;
//     offset[0]=images-bufferLen;
//     size[0]=images;
    
//     spaceDummy=DataSpace( 3, dimsf); 

//     for(int i=0; i<inputfile->outs.size();i++)
//     {
//         extendDataset(inputfile->outs[i],buffer[i],file);
//         buffer[i].resize(boost::extents[0][NX][NY]);
//     }
    
//     std::cout<<"flushing! size: "<<bufferSize/1024.0/1024<<"MB images: "<<images<<std::endl;

//     bufferLen=0;
//     bufferSize=0;
    
// }


// template<typename INTorFloat>
// void DataFile::extendDataset(std::string datasetName, const boost::multi_array<INTorFloat,3>& data_array, H5File& file)
// {  
//     DataSet dataset=file.openDataSet(datasetName);

//     dataset.extend(size);
//     DataSpace fspace = dataset.getSpace ();
//     fspace.selectHyperslab( H5S_SELECT_SET, dimsf, offset );
//     dataset.write( data_array.data(), PredType::NATIVE_INT, spaceDummy, fspace );

// }


