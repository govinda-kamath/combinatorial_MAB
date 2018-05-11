//
// Created by Govinda Kamath on 2/7/18.
//

#ifndef COMBINATORIAL_MAB_TENXREADER_H
#define COMBINATORIAL_MAB_TENXREADER_H


#include <unordered_map>
namespace tenXReader{
    std::vector<unsigned> get10xMatrixSize(std::string &fileName) ;
    void get10xDataSet(std::string &fileName, std::vector<int> &dataRead, std::vector<int> &indicesData,
                       std::vector<int> &indptrData, std::vector<unsigned> &shapeData);
    void get10xDataSet(std::string &fileName, std::vector<float> &dataRead, std::vector<int> &indicesData,
                       std::vector<int> &indptrData, std::vector<unsigned> &shapeData);
    void get10xMatrix(std::string &fileName, std::vector<std::vector<float> > &denseDataMatrix) ;
    void get10xNormalisedDenseMatrix(std::string &fileName, std::vector<std::vector<float> > &denseDataMatrix) ;
    void get10xSparseMatrix(std::string &fileName, std::vector<std::unordered_map<unsigned long, int> >
    &sparseDataMatrix);
    void get10xNormalisedSparseMatrix(std::string &fileName, std::vector<std::unordered_map<unsigned long, float> >
    &sparseNormalisedDataMatrix);
}
#endif //COMBINATORIAL_MAB_TENXREADER_H