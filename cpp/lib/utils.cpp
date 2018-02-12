//
// Created by Govinda Kamath on 2/1/18.
//
#include <vector>
#include <queue>
#include <numeric>
#include <dlib/image_io.h>
#include <dlib/image_transforms.h>
#include "utils.h"
#include <glob.h>

void utils::readImageAsVector (std::string filePath, std::vector<float> &imageVec) {

    dlib::array2d <dlib::rgb_pixel> imageRGB;
    dlib::load_image(imageRGB, filePath.c_str());
    unsigned  numColumns(imageRGB.nc()), numRows(imageRGB.nr());
    unsigned numPixels(numColumns*numRows);
    unsigned vecLength(numPixels*3);

//    std::cout << numColumns <<"\t" << numRows <<
//              "\t" << numPixels << "\t" << vecLength << std::endl;

    if (imageVec.size() != vecLength){

//        std::cout << "initialising" << std::endl;
        imageVec.clear();
        imageVec.reserve(vecLength);

        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].red);
            }
        }



        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].blue);
            }
        }

        for (unsigned i(0); i < numRows; i++){
            for (unsigned j(0); j < numColumns; j++) {
                imageVec.push_back((float) imageRGB[i][j].green);
            }
        }

    } else {

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[i * numRows + j] = (float) imageRGB[i][j].red;
            }
        }

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[numPixels + i * numRows + j] = (float) imageRGB[i][j].blue;
            }
        }

        for (unsigned i(0); i < numRows; i++) {
            for (unsigned j(0); j < numColumns; j++) {
                imageVec[2 * numPixels + i * numRows + j] = (float) imageRGB[i][j].green;
            }
        }
    }
}

void utils::getPathToFile(std::vector<std::string> & pathsToImages, const std::string directoryPath,
                   const std::string  fileSuffix){
    std::string filePrefix = "/*";
    std::string searchName = directoryPath + filePrefix + fileSuffix;
    std::cout << searchName << std::endl;
    glob_t glob_result;
    glob(searchName.c_str(), GLOB_TILDE, NULL, &glob_result);


    std::cout << "Number of files " << glob_result.gl_pathc << std::endl;

    pathsToImages.reserve(glob_result.gl_pathc);
    for (unsigned long i(0); i < glob_result.gl_pathc; i ++){
        pathsToImages.push_back(glob_result.gl_pathv[i]);
    }
    std::ofstream dataIDfile;
    dataIDfile.open ("data_id.map", std::ofstream::out | std::ofstream::app);
    for (unsigned long i(0); i< pathsToImages.size(); i++){
        dataIDfile << i << " " << pathsToImages[i] << "\n";
    }
    std::sort(pathsToImages.begin(), pathsToImages.end());
}
