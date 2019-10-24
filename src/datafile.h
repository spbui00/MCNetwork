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
//     std::string filename="out.h5";
    
//     std::vector<std::string> outs;
//     std::map<std::string,std::function<const boost::multi_array<int,2>()>> getterMap;
//     std::vector<boost::multi_array<int,3>> buffer;
//     int bufferLen=0;
//     long int bufferSize=0;
//     long int maxBufferSize=1024*1024*10;
    
    
//     unsigned int   NX = 0;          
//     unsigned int   NY = 0;
//     unsigned int images = 0;

//     hsize_t     dimsf[3];
//     hsize_t     offset[3]={0,0,0};
//     hsize_t      size[3];
//     hsize_t      chunk_dims[3];
//     hsize_t      maxdimsf[3] = {H5S_UNLIMITED, H5S_UNLIMITED,H5S_UNLIMITED};

//     DataSpace spaceDummy;
//     std::shared_ptr<ParameterStorage> parameterStorage;

    
// public:


//     ~DataFile();
//     DataFile(Lipidsystem&,std::shared_ptr<InputFile>);
//     void createFile();
//     void writeStep();

// private:
//     template<typename INTorFloat>
//     void extendDataset(std::string ,const boost::multi_array<INTorFloat,3>&,H5File&);
//     void createDataset(std::string, H5File&);
//     void createAttribute(std::string, double, H5File&);
//     void flush();
    
};

#endif // DATAFILE_H
