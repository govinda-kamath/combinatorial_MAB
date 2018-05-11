//
// Created by Govinda Kamath on 2/7/18.
//


#include<iostream>
#include<vector>
#include <cassert>
#include<H5Cpp.h>
#include <unordered_map>
#include "tenXReader.h"
#include <random>

using namespace H5;

std::vector<unsigned> tenXReader::get10xMatrixSize(std::string &fileName) {
    H5File *file = new H5File(fileName.c_str(), H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup("mm10"));

    DataSet *dataShape;
    dataShape = new DataSet(group->openDataSet("shape"));

    DataSpace dspaceShape = dataShape->getSpace();

    hsize_t rankShape;
    hsize_t dimsShape[2];

    rankShape = dspaceShape.getSimpleExtentDims(dimsShape, NULL);

    std::vector<unsigned> shapeData;
    shapeData.resize(dimsShape[0]);
    hsize_t dimsmShape[1];
    dimsmShape[0] = dimsShape[0];
    DataSpace memspaceShape(1, dimsmShape);

    dataShape->read(shapeData.data(), PredType::NATIVE_INT, memspaceShape, dspaceShape);

    return shapeData;

}

void tenXReader::get10xDataSet(std::string &fileName, std::vector<int> &dataRead, std::vector<int> &indicesData,
                               std::vector<int> &indptrData, std::vector<unsigned> &shapeData){

    H5File *file = new H5File(fileName.c_str(), H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup("mm10"));

    DataSet *dataset;
    DataSet *dataIndices;
    DataSet *dataIndptr;
    DataSet *dataShape;

    try {  // to determine if the dataset exists in the group
        dataset = new DataSet   (group->openDataSet("data"));
        dataIndices = new DataSet(group->openDataSet("indices"));
        dataIndptr = new DataSet(group->openDataSet("indptr"));
        dataShape = new DataSet(group->openDataSet("shape"));
    }
    catch (GroupIException not_found_error) {
        std::cout << " Dataset is not found." << std::endl;
    }

    DataSpace dspace = dataset->getSpace();
    DataSpace dspaceIndices = dataIndices->getSpace();
    DataSpace dspaceIndptr = dataIndptr->getSpace();
    DataSpace dspaceShape = dataShape->getSpace();


    hsize_t rankData, rankIndices, rankIndptr, rankShape;
    hsize_t dimsData[2], dimsIndices[2], dimsIndptr[2], dimsShape[2];

    rankData = dspace.getSimpleExtentDims(dimsData, NULL); // rank = 1
    rankIndices = dspaceIndices.getSimpleExtentDims(dimsIndices, NULL);
    rankIndptr = dspaceIndptr.getSimpleExtentDims(dimsIndptr, NULL);
    rankShape = dspaceShape.getSimpleExtentDims(dimsShape, NULL);



    dataRead.resize(dimsData[0]);
    shapeData.resize(dimsShape[0]);
    hsize_t dimsmShape[1];
    dimsmShape[0] = dimsShape[0];
    DataSpace memspaceShape(1, dimsmShape);

    indicesData.resize(dimsIndices[0]);
    hsize_t dimsmIndices[1];
    dimsmIndices[0] = dimsIndices[0];
    DataSpace memspaceIndices(1, dimsmIndices);


    indptrData.resize(dimsIndptr[0]);
    hsize_t dimsmIndptr[1];
    dimsmIndptr[0] = dimsIndptr[0];
    DataSpace memspaceIndptr(1, dimsmIndptr);

    dataShape->read(shapeData.data(), PredType::NATIVE_INT, memspaceShape, dspaceShape);

    dataIndices->read(indicesData.data(), PredType::NATIVE_INT, memspaceIndices, dspaceIndices);

    dataIndptr->read(indptrData.data(), PredType::NATIVE_INT, memspaceIndptr, dspaceIndptr);

    hsize_t dimsmData[1];
    dimsmData[0] = dimsData[0];
    DataSpace memspaceData(1, dimsmData);

    dataset->read(dataRead.data(), PredType::NATIVE_INT, memspaceData, dspace);
}



void tenXReader::get10xDataSet(std::string &fileName, std::vector<float> &dataRead, std::vector<int> &indicesData,
                               std::vector<int> &indptrData, std::vector<unsigned> &shapeData){

    H5File *file = new H5File(fileName.c_str(), H5F_ACC_RDWR);
    Group *group = new Group(file->openGroup("mm10"));

    DataSet *dataset;
    DataSet *dataIndices;
    DataSet *dataIndptr;
    DataSet *dataShape;

    try {  // to determine if the dataset exists in the group
        dataset = new DataSet   (group->openDataSet("data"));
        dataIndices = new DataSet(group->openDataSet("indices"));
        dataIndptr = new DataSet(group->openDataSet("indptr"));
        dataShape = new DataSet(group->openDataSet("shape"));
    }
    catch (GroupIException not_found_error) {
        std::cout << " Dataset is not found." << std::endl;
    }

    DataSpace dspace = dataset->getSpace();
    DataSpace dspaceIndices = dataIndices->getSpace();
    DataSpace dspaceIndptr = dataIndptr->getSpace();
    DataSpace dspaceShape = dataShape->getSpace();


    hsize_t rankData, rankIndices, rankIndptr, rankShape;
    hsize_t dimsData[2], dimsIndices[2], dimsIndptr[2], dimsShape[2];

    rankData = dspace.getSimpleExtentDims(dimsData, NULL); // rank = 1
    rankIndices = dspaceIndices.getSimpleExtentDims(dimsIndices, NULL);
    rankIndptr = dspaceIndptr.getSimpleExtentDims(dimsIndptr, NULL);
    rankShape = dspaceShape.getSimpleExtentDims(dimsShape, NULL);



    dataRead.resize(dimsData[0]);
    shapeData.resize(dimsShape[0]);
    hsize_t dimsmShape[1];
    dimsmShape[0] = dimsShape[0];
    DataSpace memspaceShape(1, dimsmShape);

    indicesData.resize(dimsIndices[0]);
    hsize_t dimsmIndices[1];
    dimsmIndices[0] = dimsIndices[0];
    DataSpace memspaceIndices(1, dimsmIndices);


    indptrData.resize(dimsIndptr[0]);
    hsize_t dimsmIndptr[1];
    dimsmIndptr[0] = dimsIndptr[0];
    DataSpace memspaceIndptr(1, dimsmIndptr);

    dataShape->read(shapeData.data(), PredType::NATIVE_INT, memspaceShape, dspaceShape);

    dataIndices->read(indicesData.data(), PredType::NATIVE_INT, memspaceIndices, dspaceIndices);

    dataIndptr->read(indptrData.data(), PredType::NATIVE_INT, memspaceIndptr, dspaceIndptr);

    hsize_t dimsmData[1];
    dimsmData[0] = dimsData[0];
    DataSpace memspaceData(1, dimsmData);

    dataset->read(dataRead.data(), PredType::NATIVE_FLOAT, memspaceData, dspace);
}



void tenXReader::get10xMatrix(std::string &fileName, std::vector<std::vector<float> > &denseDataMatrix) {

    std::vector<float> dataRead;
    std::vector<unsigned> shapeData;
    std::vector<int> indicesData;
    std::vector<int> indptrData;

    std::random_device rd;
    std::mt19937 g(9);
    std::normal_distribution<> dist(0, 0.0005);

    tenXReader::get10xDataSet(fileName, dataRead, indicesData, indptrData, shapeData);

    std::cout << denseDataMatrix.size() << " " << denseDataMatrix[0].size() << std::endl;
    std::cout << shapeData[0] << " " << shapeData[1] << std::endl;
    assert(denseDataMatrix.size() == shapeData[0]);
    assert(denseDataMatrix[0].size() == shapeData[1]);

    for (unsigned i(0); i < indptrData.size() - 1; i++) {
        for (unsigned j(indptrData[i]); j < indptrData[i + 1]; j++) {
            denseDataMatrix[indicesData[j]][i] = dataRead[j]+dist(g);
        }
    }
}


void tenXReader::get10xNormalisedDenseMatrix(std::string &fileName, std::vector<std::vector<float> > &denseDataMatrix) {

    std::vector<int> dataRead;
    std::vector<unsigned> shapeData;
    std::vector<int> indicesData;
    std::vector<int> indptrData;

    tenXReader::get10xDataSet(fileName, dataRead, indicesData, indptrData, shapeData);
    assert(denseDataMatrix.size() == shapeData[1]);
    assert(denseDataMatrix[0].size() == shapeData[0]);

    for (unsigned long i(0); i < indptrData.size() - 1; i++) {
        float numberOfMolecules(0.0);
        for (unsigned long j(indptrData[i]); j < indptrData[i + 1]; j++)
            numberOfMolecules += dataRead[j];

        for (unsigned long j(indptrData[i]); j < indptrData[i + 1]; j++)
            denseDataMatrix[i][indicesData[j]] = ((float)dataRead[j])/numberOfMolecules;

    }
}

void tenXReader::get10xSparseMatrix(std::string &fileName, std::vector<std::unordered_map<unsigned long, int> >
&sparseDataMatrix) {


    std::vector<int> dataRead;
    std::vector<unsigned > shapeData;
    std::vector<int> indicesData;
    std::vector<int> indptrData;

    tenXReader::get10xDataSet(fileName, dataRead, indicesData, indptrData, shapeData);

    assert(sparseDataMatrix.size() == shapeData[1]);

    for (unsigned long i(0); i < indptrData.size() - 1; i++) {
        sparseDataMatrix[i].reserve(indptrData[i + 1] - indptrData[i]);
    }

    for (unsigned long i(0); i < indptrData.size() - 1; i++) {
//        sparseDataMatrix[i].reserve(indptrData[i + 1]-indptrData[i]);
        for (unsigned long j(indptrData[i]); j < indptrData[i + 1]; j++) {
            sparseDataMatrix[i].insert( std::make_pair<unsigned long, int>(
                    (unsigned long) indicesData[j],(int)dataRead[j]));
        }
    }
}

void tenXReader::get10xNormalisedSparseMatrix(std::string &fileName, std::vector<std::unordered_map<unsigned long, float> >
&sparseNormalisedDataMatrix) {


    std::vector<int> dataRead;
    std::vector<unsigned> shapeData;
    std::vector<int> indicesData;
    std::vector<int> indptrData;

    tenXReader::get10xDataSet(fileName, dataRead, indicesData, indptrData, shapeData);

    assert(sparseNormalisedDataMatrix.size() == shapeData[1]);

    for (unsigned i(0); i < indptrData.size() - 1; i++) {
        sparseNormalisedDataMatrix[i].reserve(indptrData[i + 1] - indptrData[i]);
    }
    for (unsigned long i(0); i < indptrData.size() - 1; i++) {
        float numberOfMolecules(0.0);
        for (unsigned long j(indptrData[i]); j < indptrData[i + 1]; j++) {
            numberOfMolecules += dataRead[j];
        }
        for (unsigned long j(indptrData[i]); j < indptrData[i + 1]; j++) {
            sparseNormalisedDataMatrix[i].insert( std::make_pair<unsigned long, float>(
                    (unsigned long) indicesData[j],((float)dataRead[j])/numberOfMolecules));
        }
    }
}