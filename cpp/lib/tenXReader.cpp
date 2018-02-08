//
// Created by Govinda Kamath on 2/7/18.
//


#include<iostream>
#include<vector>
#include <cassert>
#include<H5Cpp.h>
#include <unordered_map>
#include "tenXReader.h"

using namespace H5;

void tenXReader::get10xMatrixSize(std::string fileName, std::vector<int> &sizeVect) {
    H5File *file = new H5File(fileName.c_str(), H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup("mm10"));

    DataSet *dataShape;
    dataShape = new DataSet(group->openDataSet("shape"));

    DataSpace dspaceShape = dataShape->getSpace();
    H5T_class_t type_classShape = dataShape->getTypeClass();

    hsize_t rankShape;
    hsize_t dimsShape[2];

    rankShape = dspaceShape.getSimpleExtentDims(dimsShape, NULL);
//    std::cout << "ShapeSize: " << dimsShape[0] << " rank " << rankShape << std::endl;

    std::vector<int> shapeData;
    shapeData.resize(dimsShape[0]);
    hsize_t dimsmShape[1];
    dimsmShape[0] = dimsShape[0];
    DataSpace memspaceShape(1, dimsmShape);

    dataShape->read(shapeData.data(), PredType::NATIVE_INT, memspaceShape, dspaceShape);

    std::cout << "Shape : ";
    for (int i = 0; i < shapeData.size(); i++) {
        std::cout << shapeData[i] << " ";
    }
    std::cout << std::endl;

    sizeVect[0] = shapeData[0];
    sizeVect[1] = shapeData[1];
}

void tenXReader::get10xMatrix(std::string fileName, std::vector<std::vector<int> > &denseDataMatrix) {

    H5File *file = new H5File(fileName.c_str(), H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup("mm10"));

    DataSet *dataset;
    DataSet *dataIndices;
    DataSet *dataIndptr;
    DataSet *dataShape;

    try {  // to determine if the dataset exists in the group
        dataset = new DataSet(group->openDataSet("data"));
        dataIndices = new DataSet(group->openDataSet("indices"));
        dataIndptr = new DataSet(group->openDataSet("indptr"));
        dataShape = new DataSet(group->openDataSet("shape"));
    }
    catch (GroupIException not_found_error) {
        std::cout << " Dataset is not found." << std::endl;
    }

//    std::cout << " Dataset is found." << std::endl;

    DataSpace dspace = dataset->getSpace();
    H5T_class_t type_class = dataset->getTypeClass();

    DataSpace dspaceIndices = dataIndices->getSpace();
    H5T_class_t type_classIndices = dataIndices->getTypeClass();

    DataSpace dspaceIndptr = dataIndptr->getSpace();
    H5T_class_t type_classIndptr = dataIndptr->getTypeClass();

    DataSpace dspaceShape = dataShape->getSpace();
    H5T_class_t type_classShape = dataShape->getTypeClass();


    hsize_t rankData, rankIndices, rankIndptr, rankShape;
    hsize_t dimsData[2], dimsIndices[2], dimsIndptr[2], dimsShape[2];

    rankData = dspace.getSimpleExtentDims(dimsData, NULL); // rank = 1
    rankIndices = dspaceIndices.getSimpleExtentDims(dimsIndices, NULL);
    rankIndptr = dspaceIndptr.getSimpleExtentDims(dimsIndptr, NULL);
    rankShape = dspaceShape.getSimpleExtentDims(dimsShape, NULL);

//    std::cout << "Datasize: " << dimsData[0] << " rank " << rankData << std::endl;
//    std::cout << "Indexsize: " << dimsIndices[0] << " rank " << rankIndices << std::endl;
//    std::cout << "IndPtrSize: " << dimsIndptr[0] << " rank " << rankIndptr << std::endl;
//    std::cout << "ShapeSize: " << dimsShape[0] << " rank " << rankShape << std::endl;

    std::vector<int> dataRead;
    dataRead.resize(dimsData[0]);
//    std::cout << "Vectsize: " << dataRead.size() << std::endl;

    std::vector<int> shapeData;
    shapeData.resize(dimsShape[0]);
    hsize_t dimsmShape[1];
    dimsmShape[0] = dimsShape[0];
    DataSpace memspaceShape(1, dimsmShape);

    std::vector<int> indicesData;
    indicesData.resize(dimsIndices[0]);
    hsize_t dimsmIndices[1];
    dimsmIndices[0] = dimsIndices[0];
    DataSpace memspaceIndices(1, dimsmIndices);

    std::vector<int> indptrData;
    indptrData.resize(dimsIndptr[0]);
    hsize_t dimsmIndptr[1];
    dimsmIndptr[0] = dimsIndptr[0];
    DataSpace memspaceIndptr(1, dimsmIndptr);

    dataShape->read(shapeData.data(), PredType::NATIVE_INT, memspaceShape, dspaceShape);

    dataIndices->read(indicesData.data(), PredType::NATIVE_INT, memspaceIndices, dspaceIndices);

    dataIndptr->read(indptrData.data(), PredType::NATIVE_INT, memspaceIndptr, dspaceIndptr);

//    for (int i = 0; i < shapeData.size(); i++) {
//        std::cout << shapeData[i] << " ";
//    }
//    std::cout << std::endl;

    hsize_t dimsmData[1];
    dimsmData[0] = dimsData[0];
    DataSpace memspaceData(1, dimsmData);

    dataset->read(dataRead.data(), PredType::NATIVE_INT, memspaceData, dspace);

    assert(denseDataMatrix.size() == shapeData[0]);
    assert(denseDataMatrix[0].size() == shapeData[1]);

    for (int i(0); i < indptrData.size() - 1; i++) {
        for (int j(indptrData[i]); j < indptrData[i + 1]; j++) {
            denseDataMatrix[indicesData[j]][i] = dataRead[j];
        }
    }
}

void tenXReader::get10xSparseMatrix(std::string fileName, std::vector<std::unordered_map<unsigned long, int> >
&sparseDataMatrix) {

    H5File *file = new H5File(fileName.c_str(), H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup("mm10"));

    DataSet *dataset;
    DataSet *dataIndices;
    DataSet *dataIndptr;
    DataSet *dataShape;

    try {  // to determine if the dataset exists in the group
        dataset = new DataSet(group->openDataSet("data"));
        dataIndices = new DataSet(group->openDataSet("indices"));
        dataIndptr = new DataSet(group->openDataSet("indptr"));
        dataShape = new DataSet(group->openDataSet("shape"));
    }
    catch (GroupIException not_found_error) {
        std::cout << " Dataset is not found." << std::endl;
    }

//    std::cout << " Dataset is found." << std::endl;

    DataSpace dspace = dataset->getSpace();
    H5T_class_t type_class = dataset->getTypeClass();

    DataSpace dspaceIndices = dataIndices->getSpace();
    H5T_class_t type_classIndices = dataIndices->getTypeClass();

    DataSpace dspaceIndptr = dataIndptr->getSpace();
    H5T_class_t type_classIndptr = dataIndptr->getTypeClass();

    DataSpace dspaceShape = dataShape->getSpace();
    H5T_class_t type_classShape = dataShape->getTypeClass();


    hsize_t rankData, rankIndices, rankIndptr, rankShape;
    hsize_t dimsData[2], dimsIndices[2], dimsIndptr[2], dimsShape[2];

    rankData = dspace.getSimpleExtentDims(dimsData, NULL); // rank = 1
    rankIndices = dspaceIndices.getSimpleExtentDims(dimsIndices, NULL);
    rankIndptr = dspaceIndptr.getSimpleExtentDims(dimsIndptr, NULL);
    rankShape = dspaceShape.getSimpleExtentDims(dimsShape, NULL);

//    std::cout << "Datasize: " << dimsData[0] << " rank " << rankData << std::endl;
//    std::cout << "Indexsize: " << dimsIndices[0] << " rank " << rankIndices << std::endl;
//    std::cout << "IndPtrSize: " << dimsIndptr[0] << " rank " << rankIndptr << std::endl;
//    std::cout << "ShapeSize: " << dimsShape[0] << " rank " << rankShape << std::endl;

    std::vector<int> dataRead;
    dataRead.resize(dimsData[0]);
//    std::cout << "Vectsize: " << dataRead.size() << std::endl;

    std::vector<int> shapeData;
    shapeData.resize(dimsShape[0]);
    hsize_t dimsmShape[1];
    dimsmShape[0] = dimsShape[0];
    DataSpace memspaceShape(1, dimsmShape);

    std::vector<int> indicesData;
    indicesData.resize(dimsIndices[0]);
    hsize_t dimsmIndices[1];
    dimsmIndices[0] = dimsIndices[0];
    DataSpace memspaceIndices(1, dimsmIndices);

    std::vector<int> indptrData;
    indptrData.resize(dimsIndptr[0]);
    hsize_t dimsmIndptr[1];
    dimsmIndptr[0] = dimsIndptr[0];
    DataSpace memspaceIndptr(1, dimsmIndptr);

    dataShape->read(shapeData.data(), PredType::NATIVE_INT, memspaceShape, dspaceShape);

    dataIndices->read(indicesData.data(), PredType::NATIVE_INT, memspaceIndices, dspaceIndices);

    dataIndptr->read(indptrData.data(), PredType::NATIVE_INT, memspaceIndptr, dspaceIndptr);

//    for (int i = 0; i < shapeData.size(); i++) {
//        std::cout << shapeData[i] << " ";
//    }
//    std::cout << std::endl;

    hsize_t dimsmData[1];
    dimsmData[0] = dimsData[0];
    DataSpace memspaceData(1, dimsmData);

    dataset->read(dataRead.data(), PredType::NATIVE_INT, memspaceData, dspace);


    assert(sparseDataMatrix.size() == shapeData[1]);

    for (int i(0); i < indptrData.size() - 1; i++) {
        for (int j(indptrData[i]); j < indptrData[i + 1]; j++) {
            sparseDataMatrix[i].insert( std::make_pair<unsigned long, int>(
                    (unsigned long) indicesData[j],(int)dataRead[j]));
        }
    }
}