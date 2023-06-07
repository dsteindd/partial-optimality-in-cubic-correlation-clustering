#pragma once
#ifndef DSSC_CSV_WRITER_HXX
#define DSSC_CSV_WRITER_HXX

#include <cassert>
#include <iostream>
#include <fstream>
#include <utility>


#include "dssc/math.hxx"
#include "random"
#include "andres/graph/multicut-cubic/persistency.hxx"

namespace dssc {

/// Gaussian Mixture Model
    class CSVWriter {
    public:

        // constructor for only triplet cost
        CSVWriter(
                std::string const &,
                std::vector<std::string> const &,
                std::vector<std::string> const & = {},
                std::string  = ","
        );

        void writeRow(std::vector<std::string> const & values);
        void close();

    private:
        std::ofstream out_;
        size_t numberOfColumns_;
        std::string separator_ = ",";
        std::vector<std::string> perRowFirstEntries_;

        std::string joinString(std::vector<std::string> const &);
        bool fileExists(std::string const &);
    };

    inline
    CSVWriter::CSVWriter(
            std::string const & fileName,
            std::vector<std::string> const & header,
            std::vector<std::string> const & perRowFirstEntries,
            std::string separator)
            :
            numberOfColumns_(header.size()),
            separator_(std::move(separator)),
            perRowFirstEntries_(perRowFirstEntries)
    {
        if (fileExists(fileName)){
            // append mode / do not write header
            out_ = std::ofstream (fileName, std::ios::app);
        } else {
            out_ = std::ofstream (fileName, std::ios::out);

            // write header
            if (out_.is_open()){
                std::string const headerString = joinString(header);
                out_ << headerString << std::endl;
            }
        }
    }

    inline
    void
    CSVWriter::writeRow(std::vector<std::string> const & values) {
        // assert not more values given than number of headers
        assert(numberOfColumns_ >= values.size() + perRowFirstEntries_.size());

        // fill values vector with NaNs
        std::vector<std::string> valuesToWrite(numberOfColumns_);
        for (size_t index = 0; index < perRowFirstEntries_.size(); ++index){
            valuesToWrite[index] = perRowFirstEntries_[index];
        }
        for (size_t index = perRowFirstEntries_.size(); index < perRowFirstEntries_.size() + values.size(); ++index){
            valuesToWrite[index] = values[index - perRowFirstEntries_.size()];
        }



        if (values.size() + perRowFirstEntries_.size() < numberOfColumns_){
            for (size_t index = values.size() + perRowFirstEntries_.size(); index < numberOfColumns_; ++index){
                valuesToWrite[index] = "NaN";
            }
        }

        if (out_.is_open()){
            std::string const valuesString = joinString(valuesToWrite);
            out_ << valuesString << std::endl;
        }
    }

    inline
    void CSVWriter::close() {
        out_.close();
    }

    inline
    std::string
    CSVWriter::joinString(const std::vector<std::string> & values) {
        std::string result;

        for (size_t index = 0; index < values.size() - 1; ++index){
            result += values[index] + separator_;
        }
        result += values[values.size() - 1];

        return result;
    }

    inline
    bool
    CSVWriter::fileExists(const std::string &fileName) {
        std::ifstream file(fileName.c_str());
        return file.good();
    }



}

#endif
