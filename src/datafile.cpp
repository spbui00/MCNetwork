#include "datafile.h"

std::map<std::string, int> DataFile::indexMap;

/*!
    contructs datafile class and creates new file.
    \param createNew true: old file is overwritten. \n false: dataset dimensions
   of existing file are read in, fails if file does not exist.
 */
DataFile::DataFile(std::string filename, bool createNew)
    : filename(filename)
{

    if (createNew) {
        H5File file = H5File(filename, H5F_ACC_TRUNC);
    } else {

        H5File file = H5File(filename, H5F_ACC_RDWR);

        // load datasets
        auto file_info = [](hid_t loc_id, const char* name, const H5L_info_t* linfo,
                             void* opdata) {
            indexMap[std::string(name)] = indexMap.size();
            return 0;
        };
        herr_t idx = H5Literate(file.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL,
            file_info, NULL);

        // create dataset properties in Datafile class
        for (auto const& [datasetName, datasetIndex] : indexMap) {
            DataSet dataset = file.openDataSet(datasetName);

            DataSpace dataspace = dataset.getSpace();

            int ndims = dataspace.getSimpleExtentNdims();

            hsize_t dimensions[ndims];
            dataspace.getSimpleExtentDims(dimensions, NULL);

            dims.push_back(ndims);
            dimsf.push_back(new hsize_t[ndims]);
            size.push_back(new hsize_t[ndims]);
            offset.push_back(new hsize_t[ndims]);
            chunk_dims.push_back(new hsize_t[ndims]);
            maxdimsf.push_back(new hsize_t[ndims]);

            dimsf.back()[0] = 1;
            size.back()[0] = dimensions[0];
            offset.back()[0] = dimensions[0] - 1;
            chunk_dims.back()[0] = 1;
            maxdimsf.back()[0] = H5S_UNLIMITED;
            for (size_t i = 1; i < ndims; i++) {
                dimsf.back()[i] = dimensions[i];
                size.back()[i] = dimensions[i];
                offset.back()[i] = 0;
                chunk_dims.back()[i] = dimensions[i];
                maxdimsf.back()[i] = H5S_UNLIMITED;
            }
        }
    }
}

std::unique_ptr<std::vector<double>>
DataFile::readFullDataset(std::string datasetName)
{

    H5File file = H5File(filename, H5F_ACC_RDWR);

    size_t dataSize = 1;
    for (size_t i = 0; i < dims[indexMap[datasetName]]; i++) {
        dataSize *= size[indexMap[datasetName]][i];
    }

    std::unique_ptr<std::vector<double>> data = std::make_unique<std::vector<double>>(dataSize);

    DataSet dataset = file.openDataSet(datasetName);
    dataset.read(data->data(), PredType::NATIVE_DOUBLE);

    return data;
}

std::unique_ptr<std::vector<double>>
DataFile::readDatasetSlice(std::string datasetName, int index)
{

    H5File file = H5File(filename, H5F_ACC_RDWR);

    hsize_t dataSize = 1;
    for (size_t i = 1; i < dims[indexMap[datasetName]]; i++) {
        dataSize *= size[indexMap[datasetName]][i];
    }

    std::unique_ptr<std::vector<double>> data = std::make_unique<std::vector<double>>(dataSize);

    DataSet dataset = file.openDataSet(datasetName);

    // get dataspace
    DataSpace dataspace = dataset.getSpace();
    offset[indexMap[datasetName]][0] = index; // its no problem to change the (globaly used) offset here, because
        // its resetted in addData() anyway
    dataspace.selectHyperslab(
        H5S_SELECT_SET, dimsf[indexMap[datasetName]],
        offset[indexMap[datasetName]]); // dimsf is already size of one slice

    // get memspace
    DataSpace memspace(1, &dataSize);

    dataset.read(data->data(), PredType::NATIVE_DOUBLE, memspace, dataspace);

    return data;
}

bool DataFile::checkDataSetExists(std::string datasetName)
{
    return indexMap.count(datasetName);
}

void DataFile::createDataset(std::string datasetName,
    std::vector<int> dimensions)
{

    H5File* file = new H5File(filename, H5F_ACC_RDWR);

    indexMap[datasetName] = numberOfDatasets;
    numberOfDatasets++;
    dims.push_back(dimensions.size() + 1);
    dimsf.push_back(new hsize_t[dims.back()]);
    size.push_back(new hsize_t[dims.back()]);
    offset.push_back(new hsize_t[dims.back()]);
    chunk_dims.push_back(new hsize_t[dims.back()]);
    maxdimsf.push_back(new hsize_t[dims.back()]);

    dimsf.back()[0] = 1;
    size.back()[0] = 0;
    offset.back()[0] = 0;
    chunk_dims.back()[0] = 1;
    maxdimsf.back()[0] = H5S_UNLIMITED;

    for (size_t i = 1; i < dims.back(); i++) {
        dimsf.back()[i] = dimensions[i - 1];
        size.back()[i] = dimensions[i - 1];
        offset.back()[i] = 0;
        chunk_dims.back()[i] = dimensions[i - 1];
        maxdimsf.back()[i] = H5S_UNLIMITED;
    }

    DataSpace dataspace(dims.back(), size.back(), maxdimsf.back());

    DSetCreatPropList cparms;
    cparms.setChunk(dims.back(), chunk_dims.back());
    int fill_val = 0;
    cparms.setFillValue(PredType::NATIVE_DOUBLE, &fill_val);

    DataSet dataset = file->createDataSet(datasetName, PredType::NATIVE_DOUBLE,
        dataspace, cparms);

    delete file;
}

/*!
     -- not tested --
 */
void DataFile::createAttribute(std::string attrName, double val)
{
    H5File file = H5File(filename, H5F_ACC_RDWR);

    const hsize_t dims = 1;
    DataSpace* dspace = new DataSpace(1, &dims);
    Attribute attr = file.createAttribute(attrName, PredType::NATIVE_DOUBLE, *dspace);
    delete dspace;

    attr.write(PredType::NATIVE_DOUBLE, &val);
}

/*!
    delete last entries in dataset
 */
void DataFile::shrinkDataset(std::string datasetName, int dataPointsToDelete)
{
    int index = indexMap.at(datasetName);
    size[index][0] -= dataPointsToDelete;

    H5File file;

    int waitTime = 0;
tryAgian:;
    try {
        file = H5File(filename, H5F_ACC_RDWR);
    } catch (FileIException const& e) {
        waitTime++;
        if (waitTime == 1) {
            std::cout << "------------> Unable to open file <------------"
                      << std::endl;
            std::cout << "---------------> error message: <--------------"
                      << std::endl;
            e.clearErrorStack();
            std::cout << "--------> trying again in 1 second <-----------"
                      << std::endl;
        } else {
            std::cout << "--------> already waiting for " << waitTime
                      << " seconds <----------" << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::seconds(1));
        goto tryAgian;
    }

    // get data space
    DataSet dataset = file.openDataSet(datasetName);
    dataset.extend(size[index]);
    // DataSpace fspace = dataset.getSpace ();
    // fspace.selectHyperslab( H5S_SELECT_SET, dimsf[index], offset[index]);

    // //get mem space
    // DataSpace memSpace (dims[index], dimsf[index]);

    // dataset.write( data , PredType::NATIVE_DOUBLE, memSpace, fspace );
}