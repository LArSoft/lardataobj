/// \file    raw.cxx
/// \brief   raw data utilities
/// \author  brebel@fnal.gov
/// modified by jti3@fnal.gov

#include "lardataobj/RawData/raw.h"

#include <bitset>
#include <iostream>
#include <iterator> // std::back_inserter()
#include <numeric>  // std::adjacent_difference()

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace raw {

  //----------------------------------------------------------
  void Compress(std::vector<short>& adc, raw::Compress_t compress)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);
    else if (compress == raw::kZeroHuffman) {
      unsigned int zerothreshold = 5;
      ZeroSuppression(adc, zerothreshold);
      CompressHuffman(adc);
    }
    else if (compress == raw::kZeroSuppression) {
      unsigned int zerothreshold = 5;
      ZeroSuppression(adc, zerothreshold);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }
  //----------------------------------------------------------
  void Compress(std::vector<short>& adc, raw::Compress_t compress, int& nearestneighbor)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);
    else if (compress == raw::kZeroHuffman) {
      unsigned int zerothreshold = 5;
      ZeroSuppression(adc, zerothreshold, nearestneighbor);
      CompressHuffman(adc);
    }
    else if (compress == raw::kZeroSuppression) {
      unsigned int zerothreshold = 5;
      ZeroSuppression(adc, zerothreshold, nearestneighbor);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }

  //----------------------------------------------------------
  void Compress(std::vector<short>& adc, raw::Compress_t compress, unsigned int& zerothreshold)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);

    else if (compress == raw::kZeroSuppression)
      ZeroSuppression(adc, zerothreshold);
    else if (compress == raw::kZeroHuffman) {
      ZeroSuppression(adc, zerothreshold);
      CompressHuffman(adc);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }
  //----------------------------------------------------------
  void Compress(std::vector<short>& adc,
                raw::Compress_t compress,
                unsigned int& zerothreshold,
                int& nearestneighbor)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);
    else if (compress == raw::kZeroSuppression)
      ZeroSuppression(adc, zerothreshold, nearestneighbor);
    else if (compress == raw::kZeroHuffman) {
      ZeroSuppression(adc, zerothreshold, nearestneighbor);
      CompressHuffman(adc);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }

  //----------------------------------------------------------
  void Compress(const boost::circular_buffer<std::vector<short>>& adcvec_neighbors,
                std::vector<short>& adc,
                raw::Compress_t compress,
                unsigned int& zerothreshold,
                int& nearestneighbor)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);
    else if (compress == raw::kZeroSuppression)
      ZeroSuppression(adcvec_neighbors, adc, zerothreshold, nearestneighbor);
    else if (compress == raw::kZeroHuffman) {
      ZeroSuppression(adcvec_neighbors, adc, zerothreshold, nearestneighbor);
      CompressHuffman(adc);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }

  //----------------------------------------------------------
  void Compress(std::vector<short>& adc,
                raw::Compress_t compress,
                unsigned int& zerothreshold,
                int pedestal,
                int& nearestneighbor,
                bool fADCStickyCodeFeature)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);
    else if (compress == raw::kZeroSuppression)
      ZeroSuppression(adc, zerothreshold, pedestal, nearestneighbor, fADCStickyCodeFeature);
    else if (compress == raw::kZeroHuffman) {
      ZeroSuppression(adc, zerothreshold, pedestal, nearestneighbor, fADCStickyCodeFeature);
      CompressHuffman(adc);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }

  //----------------------------------------------------------
  void Compress(const boost::circular_buffer<std::vector<short>>& adcvec_neighbors,
                std::vector<short>& adc,
                raw::Compress_t compress,
                unsigned int& zerothreshold,
                int pedestal,
                int& nearestneighbor,
                bool fADCStickyCodeFeature)
  {
    if (compress == raw::kHuffman)
      CompressHuffman(adc);
    else if (compress == raw::kZeroSuppression)
      ZeroSuppression(
        adcvec_neighbors, adc, zerothreshold, pedestal, nearestneighbor, fADCStickyCodeFeature);
    else if (compress == raw::kZeroHuffman) {
      ZeroSuppression(
        adcvec_neighbors, adc, zerothreshold, pedestal, nearestneighbor, fADCStickyCodeFeature);
      CompressHuffman(adc);
    }
    else if (compress == raw::kFibonacci) {
      CompressFibonacci(adc);
    }

    return;
  }

  //----------------------------------------------------------
  // Zero suppression function
  void ZeroSuppression(std::vector<short>& adc, unsigned int& zerothreshold)
  {
    const int adcsize = adc.size();
    const int zerothresholdsigned = zerothreshold;

    std::vector<short> zerosuppressed(adc.size());
    int maxblocks = adcsize / 2 + 1;
    std::vector<short> blockbegin(maxblocks);
    std::vector<short> blocksize(maxblocks);

    unsigned int nblocks = 0;
    unsigned int zerosuppressedsize = 0;

    int blockcheck = 0;

    for (int i = 0; i < adcsize; ++i) {
      int adc_current_value = std::abs(adc[i]);

      if (adc_current_value > zerothresholdsigned) {

        if (blockcheck == 0) {

          blockbegin[nblocks] = i;
          blocksize[nblocks] = 0;
          blockcheck = 1;
        }

        zerosuppressed[zerosuppressedsize] = adc[i];
        zerosuppressedsize++;
        blocksize[nblocks]++;

        if (i == adcsize - 1) nblocks++;
      }

      if (adc_current_value <= zerothresholdsigned && blockcheck == 1) {
        zerosuppressed[zerosuppressedsize] = adc[i];
        zerosuppressedsize++;
        blocksize[nblocks]++;
        nblocks++;
        blockcheck = 0;
      }
    }

    adc.resize(2 + nblocks + nblocks + zerosuppressedsize);

    adc[0] = adcsize; //fill first entry in adc with length of uncompressed vector
    adc[1] = nblocks;

    for (unsigned int i = 0; i < nblocks; ++i)
      adc[i + 2] = blockbegin[i];

    for (unsigned int i = 0; i < nblocks; ++i)
      adc[i + nblocks + 2] = blocksize[i];

    for (unsigned int i = 0; i < zerosuppressedsize; ++i)
      adc[i + nblocks + nblocks + 2] = zerosuppressed[i];
  }

  //----------------------------------------------------------
  // Zero suppression function which merges blocks if they are
  // within parameter nearestneighbor of each other
  void ZeroSuppression(std::vector<short>& adc, unsigned int& zerothreshold, int& nearestneighbor)
  {

    const int adcsize = adc.size();
    const int zerothresholdsigned = zerothreshold;

    std::vector<short> zerosuppressed(adcsize);
    int maxblocks = adcsize / 2 + 1;
    std::vector<short> blockbegin(maxblocks);
    std::vector<short> blocksize(maxblocks);

    int nblocks = 0;
    int zerosuppressedsize = 0;

    int blockstartcheck = 0;
    int endofblockcheck = 0;

    for (int i = 0; i < adcsize; ++i) {
      int adc_current_value = std::abs(adc[i]);

      if (blockstartcheck == 0) {
        if (adc_current_value > zerothresholdsigned) {
          if (nblocks > 0) {
            if ((i - nearestneighbor) <= (blockbegin[nblocks - 1] + blocksize[nblocks - 1] + 1)) {

              nblocks--;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
            else {
              blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
          }
          else {
            blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
            blocksize[nblocks] = i - blockbegin[nblocks] + 1;
            blockstartcheck = 1;
          }
        }
      }
      else if (blockstartcheck == 1) {
        if (adc_current_value > zerothresholdsigned) {
          blocksize[nblocks]++;
          endofblockcheck = 0;
        }
        else {
          if (endofblockcheck < nearestneighbor) {
            endofblockcheck++;
            blocksize[nblocks]++;
          }
          //block has ended
          else if (i + 2 < adcsize) { //check if end of adc vector is near
            if (std::abs(adc[i + 1]) <= zerothresholdsigned &&
                std::abs(adc[i + 2]) <= zerothresholdsigned) {
              endofblockcheck = 0;
              blockstartcheck = 0;
              nblocks++;
            }
          }

        } // end else
      }   // end if blockstartcheck == 1
    }     // end loop over adc size

    if (blockstartcheck == 1) { // we reached the end of the adc vector with the block still going
      ++nblocks;
    }

    for (int i = 0; i < nblocks; ++i)
      zerosuppressedsize += blocksize[i];

    adc.resize(2 + nblocks + nblocks + zerosuppressedsize);
    zerosuppressed.resize(2 + nblocks + nblocks + zerosuppressedsize);

    int zerosuppressedcount = 0;
    for (int i = 0; i < nblocks; ++i) {
      //zerosuppressedsize += blocksize[i];
      for (int j = 0; j < blocksize[i]; ++j) {
        zerosuppressed[zerosuppressedcount] = adc[blockbegin[i] + j];
        zerosuppressedcount++;
      }
    }

    adc[0] = adcsize; //fill first entry in adc with length of uncompressed vector
    adc[1] = nblocks;
    for (int i = 0; i < nblocks; ++i) {
      adc[i + 2] = blockbegin[i];
      adc[i + nblocks + 2] = blocksize[i];
    }

    for (int i = 0; i < zerosuppressedsize; ++i)
      adc[i + nblocks + nblocks + 2] = zerosuppressed[i];

    // for(int i = 0; i < 2 + 2*nblocks + zerosuppressedsize; ++i)
    //   std::cout << adc[i] << std::endl;
    //adc.resize(2+nblocks+nblocks+zerosuppressedsize);
  }

  //----------------------------------------------------------
  // Zero suppression function which merges blocks if they are
  // within parameter nearestneighbor of each other
  // after subtracting pedestal value
  void ZeroSuppression(std::vector<short>& adc,
                       unsigned int& zerothreshold,
                       int pedestal,
                       int& nearestneighbor,
                       bool fADCStickyCodeFeature)
  {

    const int adcsize = adc.size();
    const int zerothresholdsigned = zerothreshold;

    std::vector<short> zerosuppressed(adcsize);
    int maxblocks = adcsize / 2 + 1;
    std::vector<short> blockbegin(maxblocks);
    std::vector<short> blocksize(maxblocks);

    int nblocks = 0;
    int zerosuppressedsize = 0;

    int blockstartcheck = 0;
    int endofblockcheck = 0;

    for (int i = 0; i < adcsize; ++i) {
      int adc_current_value = ADCStickyCodeCheck(adc[i], pedestal, fADCStickyCodeFeature);

      if (blockstartcheck == 0) {
        if (adc_current_value > zerothresholdsigned) {
          if (nblocks > 0) {
            if (i - nearestneighbor <= blockbegin[nblocks - 1] + blocksize[nblocks - 1] + 1) {
              nblocks--;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
            else {
              blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
          }
          else {
            blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
            blocksize[nblocks] = i - blockbegin[nblocks] + 1;
            blockstartcheck = 1;
          }
        }
      }
      else if (blockstartcheck == 1) {
        if (adc_current_value > zerothresholdsigned) {
          blocksize[nblocks]++;
          endofblockcheck = 0;
        }
        else {
          if (endofblockcheck < nearestneighbor) {
            endofblockcheck++;
            blocksize[nblocks]++;
          }
          //block has ended
          else if (i + 2 < adcsize) { //check if end of adc vector is near
            if (ADCStickyCodeCheck(adc[i + 1], pedestal, fADCStickyCodeFeature) <=
                  zerothresholdsigned &&
                ADCStickyCodeCheck(adc[i + 2], pedestal, fADCStickyCodeFeature) <=
                  zerothresholdsigned) {
              endofblockcheck = 0;
              blockstartcheck = 0;
              nblocks++;
            }
          }
        } // end else
      }   // end if blockstartcheck == 1
    }     // end loop over adc size

    if (blockstartcheck == 1) { // we reached the end of the adc vector with the block still going
      ++nblocks;
    }

    for (int i = 0; i < nblocks; ++i)
      zerosuppressedsize += blocksize[i];

    adc.resize(2 + nblocks + nblocks + zerosuppressedsize);
    zerosuppressed.resize(2 + nblocks + nblocks + zerosuppressedsize);

    int zerosuppressedcount = 0;
    for (int i = 0; i < nblocks; ++i) {
      //zerosuppressedsize += blocksize[i];
      for (int j = 0; j < blocksize[i]; ++j) {
        zerosuppressed[zerosuppressedcount] = adc[blockbegin[i] + j];
        zerosuppressedcount++;
      }
    }

    adc[0] = adcsize; //fill first entry in adc with length of uncompressed vector
    adc[1] = nblocks;
    for (int i = 0; i < nblocks; ++i) {
      adc[i + 2] = blockbegin[i];
      adc[i + nblocks + 2] = blocksize[i];
    }

    for (int i = 0; i < zerosuppressedsize; ++i)
      adc[i + nblocks + nblocks + 2] = zerosuppressed[i];

    // for(int i = 0; i < 2 + 2*nblocks + zerosuppressedsize; ++i)
    //   std::cout << adc[i] << std::endl;
    //adc.resize(2+nblocks+nblocks+zerosuppressedsize);
  }

  //----------------------------------------------------------
  // Zero suppression function which merges blocks if they are
  // within parameter nearest neighbor of each other and makes
  // blocks if neighboring wires have nonzero blocks there
  void ZeroSuppression(const boost::circular_buffer<std::vector<short>>& adcvec_neighbors,
                       std::vector<short>& adc,
                       unsigned int& zerothreshold,
                       int& nearestneighbor)
  {

    const int adcsize = adc.size();
    const int zerothresholdsigned = zerothreshold;

    std::vector<short> zerosuppressed(adcsize);
    const int maxblocks = adcsize / 2 + 1;
    std::vector<short> blockbegin(maxblocks);
    std::vector<short> blocksize(maxblocks);

    int nblocks = 0;
    int zerosuppressedsize = 0;

    int blockstartcheck = 0;
    int endofblockcheck = 0;

    for (int i = 0; i < adcsize; ++i) {

      //find maximum adc value among all neighboring channels within the nearest neighbor channel distance

      int adc_current_value = 0;

      for (boost::circular_buffer<std::vector<short>>::const_iterator adcveciter =
             adcvec_neighbors.begin();
           adcveciter != adcvec_neighbors.end();
           ++adcveciter) {
        const std::vector<short>& adcvec_current = *adcveciter;
        const int adcvec_current_single = std::abs(adcvec_current[i]);

        if (adc_current_value < adcvec_current_single) {
          adc_current_value = adcvec_current_single;
        }
      }
      if (blockstartcheck == 0) {
        if (adc_current_value > zerothresholdsigned) {
          if (nblocks > 0) {
            if (i - nearestneighbor <= blockbegin[nblocks - 1] + blocksize[nblocks - 1] + 1) {
              nblocks--;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
            else {
              blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
          }
          else {
            blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
            blocksize[nblocks] = i - blockbegin[nblocks] + 1;
            blockstartcheck = 1;
          }
        }
      }
      else if (blockstartcheck == 1) {
        if (adc_current_value > zerothresholdsigned) {
          blocksize[nblocks]++;
          endofblockcheck = 0;
        }
        else {
          if (endofblockcheck < nearestneighbor) {
            endofblockcheck++;
            blocksize[nblocks]++;
          }
          //block has ended
          else if (i + 2 < adcsize) { //check if end of adc vector is near
            if (std::abs(adc[i + 1]) <= zerothresholdsigned &&
                std::abs(adc[i + 2]) <= zerothresholdsigned) {
              endofblockcheck = 0;
              blockstartcheck = 0;
              nblocks++;
            }
          }

        } // end else
      }   // end if blockstartcheck == 1
    }     // end loop over adc size

    if (blockstartcheck == 1) { // we reached the end of the adc vector with the block still going
      ++nblocks;
    }

    for (int i = 0; i < nblocks; ++i)
      zerosuppressedsize += blocksize[i];

    adc.resize(2 + nblocks + nblocks + zerosuppressedsize);
    zerosuppressed.resize(2 + nblocks + nblocks + zerosuppressedsize);

    int zerosuppressedcount = 0;
    for (int i = 0; i < nblocks; ++i) {
      //zerosuppressedsize += blocksize[i];
      for (int j = 0; j < blocksize[i]; ++j) {
        zerosuppressed[zerosuppressedcount] = adc[blockbegin[i] + j];
        zerosuppressedcount++;
      }
    }

    adc[0] = adcsize; //fill first entry in adc with length of uncompressed vector
    adc[1] = nblocks;
    for (int i = 0; i < nblocks; ++i) {
      adc[i + 2] = blockbegin[i];
      adc[i + nblocks + 2] = blocksize[i];
    }

    for (int i = 0; i < zerosuppressedsize; ++i)
      adc[i + nblocks + nblocks + 2] = zerosuppressed[i];

    // for(int i = 0; i < 2 + 2*nblocks + zerosuppressedsize; ++i)
    //   std::cout << adc[i] << std::endl;
    //adc.resize(2+nblocks+nblocks+zerosuppressedsize);
  }

  //----------------------------------------------------------
  // Zero suppression function which merges blocks if they are
  // within parameter nearest neighbor of each other and makes
  // blocks if neighboring wires have nonzero blocks there
  // after subtracting pedestal values
  void ZeroSuppression(const boost::circular_buffer<std::vector<short>>& adcvec_neighbors,
                       std::vector<short>& adc,
                       unsigned int& zerothreshold,
                       int pedestal,
                       int& nearestneighbor,
                       bool fADCStickyCodeFeature)
  {

    const int adcsize = adc.size();
    const int zerothresholdsigned = zerothreshold;

    std::vector<short> zerosuppressed(adcsize);
    const int maxblocks = adcsize / 2 + 1;
    std::vector<short> blockbegin(maxblocks);
    std::vector<short> blocksize(maxblocks);

    int nblocks = 0;
    int zerosuppressedsize = 0;

    int blockstartcheck = 0;
    int endofblockcheck = 0;

    for (int i = 0; i < adcsize; ++i) {

      //find maximum adc value among all neighboring channels within the nearest neighbor channel distance

      int adc_current_value = ADCStickyCodeCheck(adc[i], pedestal, fADCStickyCodeFeature);

      for (boost::circular_buffer<std::vector<short>>::const_iterator adcveciter =
             adcvec_neighbors.begin();
           adcveciter != adcvec_neighbors.end();
           ++adcveciter) {
        const std::vector<short>& adcvec_current = *adcveciter;
        const int adcvec_current_single = std::abs(adcvec_current[i] - pedestal);

        if (adc_current_value < adcvec_current_single) {
          adc_current_value = adcvec_current_single;
        }
      }
      if (blockstartcheck == 0) {
        if (adc_current_value > zerothresholdsigned) {
          if (nblocks > 0) {
            if (i - nearestneighbor <= blockbegin[nblocks - 1] + blocksize[nblocks - 1] + 1) {
              nblocks--;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
            else {
              blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
              blocksize[nblocks] = i - blockbegin[nblocks] + 1;
              blockstartcheck = 1;
            }
          }
          else {
            blockbegin[nblocks] = (i - nearestneighbor > 0) ? i - nearestneighbor : 0;
            blocksize[nblocks] = i - blockbegin[nblocks] + 1;
            blockstartcheck = 1;
          }
        }
      }
      else if (blockstartcheck == 1) {
        if (adc_current_value > zerothresholdsigned) {
          blocksize[nblocks]++;
          endofblockcheck = 0;
        }
        else {
          if (endofblockcheck < nearestneighbor) {
            endofblockcheck++;
            blocksize[nblocks]++;
          }
          //block has ended

          else if (i + 2 < adcsize) { //check if end of adc vector is near
            if (ADCStickyCodeCheck(adc[i + 1], pedestal, fADCStickyCodeFeature) <=
                  zerothresholdsigned &&
                ADCStickyCodeCheck(adc[i + 2], pedestal, fADCStickyCodeFeature) <=
                  zerothresholdsigned) {
              endofblockcheck = 0;
              blockstartcheck = 0;
              nblocks++;
            }
          }

        } // end else
      }   // end if blockstartcheck == 1
    }     // end loop over adc size

    if (blockstartcheck == 1) { // we reached the end of the adc vector with the block still going
      ++nblocks;
    }

    for (int i = 0; i < nblocks; ++i)
      zerosuppressedsize += blocksize[i];

    adc.resize(2 + nblocks + nblocks + zerosuppressedsize);
    zerosuppressed.resize(2 + nblocks + nblocks + zerosuppressedsize);

    int zerosuppressedcount = 0;
    for (int i = 0; i < nblocks; ++i) {
      //zerosuppressedsize += blocksize[i];
      for (int j = 0; j < blocksize[i]; ++j) {
        zerosuppressed[zerosuppressedcount] = adc[blockbegin[i] + j];
        zerosuppressedcount++;
      }
    }

    adc[0] = adcsize; //fill first entry in adc with length of uncompressed vector
    adc[1] = nblocks;
    for (int i = 0; i < nblocks; ++i) {
      adc[i + 2] = blockbegin[i];
      adc[i + nblocks + 2] = blocksize[i];
    }

    for (int i = 0; i < zerosuppressedsize; ++i)
      adc[i + nblocks + nblocks + 2] = zerosuppressed[i];

    // for(int i = 0; i < 2 + 2*nblocks + zerosuppressedsize; ++i)
    //   std::cout << adc[i] << std::endl;
    //adc.resize(2+nblocks+nblocks+zerosuppressedsize);
  }

  //----------------------------------------------------------
  // Reverse zero suppression function
  void ZeroUnsuppression(const std::vector<short>& adc, std::vector<short>& uncompressed)
  {
    const int lengthofadc = adc[0];
    const int nblocks = adc[1];

    uncompressed.resize(lengthofadc);
    for (int i = 0; i < lengthofadc; ++i) {
      uncompressed[i] = 0;
    }

    int zerosuppressedindex = nblocks * 2 + 2;

    for (int i = 0; i < nblocks; ++i) { //loop over each nonzero block of the compressed vector

      for (int j = 0; j < adc[2 + nblocks + i]; ++j) { //loop over each block size

        //set uncompressed value
        uncompressed[adc[2 + i] + j] = adc[zerosuppressedindex];
        zerosuppressedindex++;
      }
    }

    return;
  }

  //----------------------------------------------------------
  // Reverse zero suppression function with pedestal re-addition
  void ZeroUnsuppression(const std::vector<short>& adc,
                         std::vector<short>& uncompressed,
                         int pedestal)
  {
    const int lengthofadc = adc[0];
    const int nblocks = adc[1];

    uncompressed.resize(lengthofadc);
    for (int i = 0; i < lengthofadc; ++i) {
      uncompressed[i] = pedestal;
    }

    int zerosuppressedindex = nblocks * 2 + 2;

    for (int i = 0; i < nblocks; ++i) { //loop over each nonzero block of the compressed vector

      for (int j = 0; j < adc[2 + nblocks + i]; ++j) { //loop over each block size

        //set uncompressed value
        uncompressed[adc[2 + i] + j] = adc[zerosuppressedindex];
        zerosuppressedindex++;
      }
    }

    return;
  }

  //----------------------------------------------------------
  // if the compression type is kNone, copy the adc vector into the uncompressed vector
  void Uncompress(const std::vector<short>& adc,
                  std::vector<short>& uncompressed,
                  raw::Compress_t compress)
  {
    if (compress == raw::kHuffman)
      UncompressHuffman(adc, uncompressed);
    else if (compress == raw::kZeroSuppression) {
      ZeroUnsuppression(adc, uncompressed);
    }
    else if (compress == raw::kZeroHuffman) {
      std::vector<short> tmp(2 * adc[0]);
      UncompressHuffman(adc, tmp);
      ZeroUnsuppression(tmp, uncompressed);
    }
    else if (compress == raw::kNone) {
      for (unsigned int i = 0; i < adc.size(); ++i)
        uncompressed[i] = adc[i];
    }
    else if (compress == raw::kFibonacci) {
      UncompressFibonacci(adc, uncompressed);
    }
    else {
      throw cet::exception("raw") << "raw::Uncompress() does not support compression #"
                                  << ((int)compress);
    }
    return;
  }

  //----------------------------------------------------------
  // if the compression type is kNone, copy the adc vector into the uncompressed vector
  void Uncompress(const std::vector<short>& adc,
                  std::vector<short>& uncompressed,
                  int pedestal,
                  raw::Compress_t compress)
  {
    if (compress == raw::kHuffman)
      UncompressHuffman(adc, uncompressed);
    else if (compress == raw::kZeroSuppression) {
      ZeroUnsuppression(adc, uncompressed, pedestal);
    }
    else if (compress == raw::kZeroHuffman) {
      std::vector<short> tmp(2 * adc[0]);
      UncompressHuffman(adc, tmp);
      ZeroUnsuppression(tmp, uncompressed, pedestal);
    }
    else if (compress == raw::kNone) {
      for (unsigned int i = 0; i < adc.size(); ++i)
        uncompressed[i] = adc[i];
    }
    else if (compress == raw::kFibonacci) {
      UncompressFibonacci(adc, uncompressed);
    }
    else {
      throw cet::exception("raw") << "raw::Uncompress() does not support compression #"
                                  << ((int)compress);
    }
    return;
  }

  // the current Huffman Coding scheme used by uBooNE is
  // based on differences between adc values in adjacent time bins
  // the code is
  // no change for 4 ticks --> 1
  // no change for 1 tick  --> 01
  // +1 change             --> 001
  // -1 change             --> 0001
  // +2 change             --> 00001
  // -2 change             --> 000001
  // +3 change             --> 0000001
  // -3 change             --> 00000001
  // abs(change) > 3       --> write actual value to short
  // use 15th bit to set whether a block is encoded or raw value
  // 1 --> Huffman coded, 0 --> raw
  // pad out the lowest bits in a word with 0's
  void CompressHuffman(std::vector<short>& adc)
  {
    std::vector<short> const orig_adc(std::move(adc));

    // diffs contains the difference between an element of adc and the previous
    // one; the first entry is never used.
    std::vector<short> diffs;
    diffs.reserve(orig_adc.size());
    std::adjacent_difference(orig_adc.begin(), orig_adc.end(), std::back_inserter(diffs));

    // prepare adc for the new data; we kind-of-expect the size,
    // so we pre-allocate it; we might want to shrink-to-fit at the end
    adc.clear();
    adc.reserve(orig_adc.size());
    // now loop over the diffs and do the Huffman encoding
    adc.push_back(orig_adc.front());
    unsigned int curb = 15U;

    std::bitset<16> bset;
    bset.set(15);

    for (size_t i = 1U; i < diffs.size(); ++i) {

      switch (diffs[i]) {
      // if the difference is 0, check to see what the next 3 differences are
      case 0: {
        if (i < diffs.size() - 3) {
          // if next 3 are also 0, set the next bit to be 1
          if (diffs[i + 1] == 0 && diffs[i + 2] == 0 && diffs[i + 3] == 0) {
            if (curb > 0) {
              --curb;
              bset.set(curb);
              i += 3;
              continue;
            }
            else {
              adc.push_back(bset.to_ulong());

              // reset the bitset to be ready for the next word
              bset.reset();
              bset.set(15);
              bset.set(14); // account for the fact that this is a zero diff
              curb = 14;
              i += 3;
              continue;
            } // end if curb is not big enough to put current difference in bset
          }   // end if next 3 are also zero
          else {
            // 0 diff is encoded as 01, so move the current bit one to the right
            if (curb > 1) {
              curb -= 2;
              bset.set(curb);
              continue;
            } // end if the current bit is large enough to set this one
            else {
              adc.push_back(bset.to_ulong());
              // reset the bitset to be ready for the next word
              bset.reset();
              bset.set(15);
              bset.set(13); // account for the fact that this is a zero diff
              curb = 13;
              continue;
            } // end if curb is not big enough to put current difference in bset
          }   // end if next 3 are not also 0
        }     // end if able to check next 3
        else {
          // 0 diff is encoded as 01, so move the current bit one to the right
          if (curb > 1) {
            curb -= 2;
            bset.set(curb);
            continue;
          } // end if the current bit is large enough to set this one
          else {
            adc.push_back(bset.to_ulong());
            // reset the bitset to be ready for the next word
            bset.reset();
            bset.set(15);
            bset.set(13); // account for the fact that this is a zero diff
            curb = 13;
            continue;
          } // end if curb is not big enough to put current difference in bset
        }   // end if not able to check the next 3
        break;
      } // end if current difference is zero
      case 1: {
        if (curb > 2) {
          curb -= 3;
          bset.set(curb);
        }
        else {
          adc.push_back(bset.to_ulong());
          // reset the bitset to be ready for the next word
          bset.reset();
          bset.set(15);
          bset.set(12); // account for the fact that this is a +1 diff
          curb = 12;
        } // end if curb is not big enough to put current difference in bset
        break;
      } // end if difference = 1
      case -1: {
        if (curb > 3) {
          curb -= 4;
          bset.set(curb);
        }
        else {
          adc.push_back(bset.to_ulong());
          // reset the bitset to be ready for the next word
          bset.reset();
          bset.set(15);
          bset.set(11); // account for the fact that this is a -1 diff
          curb = 11;
        } // end if curb is not big enough to put current difference in bset
        break;
      } // end if difference = -1
      case 2: {
        if (curb > 4) {
          curb -= 5;
          bset.set(curb);
        }
        else {
          adc.push_back(bset.to_ulong());
          // reset the bitset to be ready for the next word
          bset.reset();
          bset.set(15);
          bset.set(10); // account for the fact that this is a +2 diff
          curb = 10;
        } // end if curb is not big enough to put current difference in bset
        break;
      } // end if difference = 2
      case -2: {
        if (curb > 5) {
          curb -= 6;
          bset.set(curb);
        }
        else {
          adc.push_back(bset.to_ulong());
          // reset the bitset to be ready for the next word
          bset.reset();
          bset.set(15);
          bset.set(9); // account for the fact that this is a -2 diff
          curb = 9;
        } // end if curb is not big enough to put current difference in bset
        break;
      } // end if difference = -2
      case 3: {
        if (curb > 6) {
          curb -= 7;
          bset.set(curb);
        }
        else {
          adc.push_back(bset.to_ulong());
          // reset the bitset to be ready for the next word
          bset.reset();
          bset.set(15);
          bset.set(8); // account for the fact that this is a +3 diff
          curb = 8;
        } // end if curb is not big enough to put current difference in bset
        break;
      } // end if difference = 3
      case -3: {
        if (curb > 7) {
          curb -= 8;
          bset.set(curb);
        }
        else {
          adc.push_back(bset.to_ulong());
          // reset the bitset to be ready for the next word
          bset.reset();
          bset.set(15);
          bset.set(7); // account for the fact that this is a -3 diff
          curb = 7;
        } // end if curb is not big enough to put current difference in bset
        break;
      } // end if difference = -3
      default: {
        // if the difference is too large that we have to put the entire adc value in:
        // put the current value into the adc vec unless the current bit is 15, then there
        // were multiple large difference values in a row
        if (curb != 15) { adc.push_back(bset.to_ulong()); }

        bset.reset();
        bset.set(15);
        curb = 15;

        // put the current adcvalue in adc, with its bit 15 set to 0
        if (orig_adc[i] > 0)
          adc.push_back(orig_adc[i]);
        else {
          std::bitset<16> tbit(-orig_adc[i]);
          tbit.set(14);
          adc.push_back(tbit.to_ulong());
        }
        break;
      } // if |difference| > 3
      } // switch diff[i]
    }   // end loop over differences

    //write out the last bitset
    adc.push_back(bset.to_ulong());

    // this would reduce global memory usage,
    // at the cost of a new allocation and copy
    //  adc.shrink_to_fit();

  } // CompressHuffman()
  //--------------------------------------------------------
  // need to decrement the bit you are looking at to determine the deltas as that is how
  // the bits are set
  void UncompressHuffman(const std::vector<short>& adc, std::vector<short>& uncompressed)
  {

    //the first entry in adc is a data value by construction
    uncompressed[0] = adc[0];

    unsigned int curu = 1;
    short curADC = uncompressed[0];

    // loop over the entries in adc and uncompress them according to the
    // encoding scheme above the CompressHuffman method
    for (unsigned int i = 1; i < adc.size() && curu < uncompressed.size(); ++i) {

      std::bitset<16> bset(adc[i]);

      int numu = 0;

      //check the 15 bit to see if this entry is a full data value or not
      if (!bset.test(15)) {
        curADC = adc[i];
        if (bset.test(14)) {
          bset.set(14, false);
          curADC = -1 * bset.to_ulong();
        }
        uncompressed[curu] = curADC;

        ++curu;
      }
      else {

        int b = 14;
        int lowestb = 0;

        // ignore any padding with zeros in the lower order bits
        while (!bset.test(lowestb) && lowestb < 15)
          ++lowestb;

        if (lowestb > 14) {
          mf::LogWarning("raw.cxx")
            << "encoded entry has no set bits!!! " << i << " "
            << bset.to_string<char, std::char_traits<char>, std::allocator<char>>();
          continue;
        }

        while (b >= lowestb) {

          // count the zeros between the current bit and the next on bit
          int zerocnt = 0;
          while (!bset.test(b - zerocnt) && b - zerocnt > lowestb)
            ++zerocnt;

          b -= zerocnt;

          if (zerocnt == 0) {
            for (int s = 0; s < 4; ++s) {
              uncompressed[curu] = curADC;
              ++curu;
              ++numu;
              if (curu > uncompressed.size() - 1) break;
            }
            --b;
          }
          else if (zerocnt == 1) {
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }
          else if (zerocnt == 2) {
            curADC += 1;
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }
          else if (zerocnt == 3) {
            curADC -= 1;
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }
          else if (zerocnt == 4) {
            curADC += 2;
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }
          else if (zerocnt == 5) {
            curADC -= 2;
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }
          else if (zerocnt == 6) {
            curADC += 3;
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }
          else if (zerocnt == 7) {
            curADC -= 3;
            uncompressed[curu] = curADC;
            ++curu;
            ++numu;
            --b;
          }

          if (curu > uncompressed.size() - 1) break;

        } // end loop over bits

        if (curu > uncompressed.size() - 1) break;

      } // end if this entry in the vector is encoded

    } // end loop over entries in adc

    return;
  }

  //--------------------------------------------------------
  // need to decrement the bit you are looking at to determine the deltas as that is how
  // the bits are set
  int ADCStickyCodeCheck(const short adc_value, const int pedestal, bool fADCStickyCodeFeature)
  {

    int adc_return_value = std::abs(adc_value - pedestal);

    if (!fADCStickyCodeFeature) { return adc_return_value; }
    // if DUNE 35t ADC sticky code feature is enabled in simulation, skip over ADC codes with LSBs of 0x00 or 0x3f
    unsigned int sixlsbs = adc_value & onemask;

    if ((sixlsbs == onemask || sixlsbs == 0) && std::abs(adc_value - pedestal) < 64) {
      adc_return_value =
        0; //set current adc value to zero if its LSBs are at sticky values and if it is within one MSB cell (64 ADC counts) of the pedestal value
    }
    return adc_return_value;
  }

  inline void add_to_sequence_terminate(std::vector<bool> array, std::vector<bool>& cmp)
  {
    cmp.insert(cmp.end(), array.begin(), array.end());
    cmp.push_back(1);
  }

  void CompressFibonacci(std::vector<short>& wf,
                         std::function<void(int, std::vector<std::vector<bool>>&)> add_to_table)
  {

    std::vector<std::vector<bool>> table;
    add_to_table(100, table);

    std::vector<short> comp_short;
    // First numbers are not encoded (size and baseline)
    assert(not empty(wf));

    unsigned int size_short = sizeof(wf[0]) * 8 - 1; // assuming we don't encode over the sign bits
    unsigned int max_2_short =
      (1 << (size_short *
             2)); //this number is then: 1,073,741,824 (i.e. almost 9min of data at 2MHz)
    size_t wf_size = wf.size();

    if (wf_size > max_2_short) {
      throw cet::exception("raw")
        << "WOW! You are trying to compress a " << wf_size << " long waveform"
        << " in a vector of shorts.\nUnfortunately, the encoded waveform needs to store the size"
        << " of the wavefom (in its first number) and " << wf_size
        << " is bigger than the maximum\n"
        << " number you can encode with 2 shorts (" << max_2_short << ").\n"
        << " Bailing out disgracefully to avoid massive trouble.\n";
    }

    short high = (wf_size >> (sizeof(short) * 8 - 1));
    short low = wf_size % ((std::numeric_limits<short>::max() + 1));
    comp_short.push_back(high);
    comp_short.push_back(low);
    comp_short.push_back(*wf.begin());

    // The format we are working with
    std::vector<bool> cmp;
    cmp.reserve(wf.size() * sizeof(wf[0]) * 8);

    // The input is changed to be the difference between ticks
    std::vector<short> diff;
    diff.reserve(wf.size());
    // Fibonacci number are positive, so need a way to get negative number
    // We use the ZigZag approach (even numbers are positive, odd are negative)
    std::adjacent_difference(
      wf.begin(), wf.end(), std::back_inserter(diff), [](const short& it1, const short& it2) {
        short d = it1 - it2;
        if (d > 0)
          d = 2 * d;
        else
          d = -2 * d + 1;
        return d;
      });

    // Start from 1 to avoid the first number
    for (size_t iSample = 1; iSample < diff.size(); ++iSample) {

      short d = diff[iSample];
      if ((unsigned)d < table.size()) { add_to_sequence_terminate(table[d], cmp); }
      else { // catch if the table is too small...
        int end = d + 1;
        // use the user provided function to fill the table
        add_to_table(end, table);
        // and add again to the sequence
        add_to_sequence_terminate(table[d], cmp);
      }
    }

    // Now convert all this to a vector of short and it's just another day in paradise for larsoft
    size_t n_vector = cmp.size();

    // Create a bitset of the size of the short to simplify things
    std::bitset<8 * sizeof(short)> this_encoded;
    // Set the bitset to match the stream
    size_t bit_counter = 0;

    for (size_t it = 0; it < n_vector; ++it) {

      if (bit_counter >= 8 * sizeof(short)) {
        short comp_s = (short)this_encoded.to_ulong();
        comp_short.push_back(comp_s);
        bit_counter = 0;
        this_encoded.reset();
      }

      if (cmp[it]) this_encoded.set(bit_counter);

      bit_counter++;
    }

    // Deal with the last part
    short comp_s = (short)this_encoded.to_ulong();
    comp_short.push_back(comp_s);

    wf = comp_short;

    return;
  }

  void UncompressFibonacci(const std::vector<short>& adc,
                           std::vector<short>& uncompressed,
                           std::function<int(std::vector<bool>&)> decode_table_chunk)
  {

    // First compressed sample is the size
    size_t n_samples = (adc[0] << (sizeof(short) * 8 - 1)) + adc[1];
    // The second compressed sample is the first uncompressed sample
    uncompressed.push_back(adc[2]);

    // The thing that we want to decode (rather than jumbled short vector)
    std::vector<bool> comp;

    for (size_t i = 3; i < adc.size(); ++i) {

      std::bitset<8 * sizeof(short)> this_encoded(adc[i]);
      for (size_t i2 = 0; i2 < this_encoded.size(); ++i2) {
        comp.push_back(this_encoded[i2]);
      }
    }

    // The bit which has to be decoded ("chunk")
    std::vector<bool> current_number;
    for (size_t it = 0; it < comp.size(); ++it) {
      current_number.push_back(comp[it]);
      // If we have at least 2 numbers in the current chunk
      // and if the last 2 number are one
      // and if we are not at the end
      if ((current_number.size() >= 2 && it >= 1 && comp[it - 1] == 1 && comp[it] == 1) ||
          it == comp.size() - 1) {
        current_number.pop_back();
        short zigzag_number = decode_table_chunk(current_number);

        short decoded = 0;
        if (zigzag_number % 2 == 0)
          decoded = zigzag_number / 2;
        else
          decoded = -(zigzag_number - 1) / 2;
        short baseline = uncompressed.back();

        uncompressed.push_back(baseline + decoded);

        current_number.clear();
        if (uncompressed.size() == n_samples) break;
      }
    }
    return;
  }

  void fibonacci_encode_table(int end, std::vector<std::vector<bool>>& table)
  {
    if ((int)table.size() > end) return;

    std::map<int, int> fibn; // no need for this to be static
    fibn[0] = 0;
    fibn[1] = 1;
    fibn[2] = 1;

    // if the table was empty, fill it
    if (empty(table)) {
      table.push_back(
        std::vector<
          bool>()); // we need this one so that the index of the vector is the same as the number we want to encode
      table.push_back({true});
      table.push_back({false, true});
    }

    if ((int)table.size() > end) return;

    // start at the third fibonacci number
    for (int i = 2; fibn.rbegin()->second < end; ++i) {
      fibn[i] = fibn[i - 1] + fibn[i - 2];
    }

    for (int i = table.size(); i <= end; ++i) {
      int current_number = i;
      std::vector<bool> seq; // the final sequence of numbers
      while (current_number > 0) {
        // find the biggest fibonacci number that is smaller than the current number
        for (auto fib = fibn.rbegin(); fib != fibn.rend(); ++fib) {
          if (fib->second <= current_number) {
            // if this is the first number, create the vector
            if (!seq.size()) seq = std::vector<bool>(fib->first - 1, false);
            // fill the sequence
            seq[fib->first - 2] = true;
            // subtract the current number
            current_number -= fib->second;
            break;
          }
        }
      }

      table.push_back(seq);
    }
  }

  short fibonacci_decode(std::vector<bool>& chunk)
  {
    std::vector<int> FibNumbers = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987};

    if (chunk.size() > FibNumbers.size()) {
      for (int i = FibNumbers.size(); i <= (int)chunk.size(); ++i) {
        FibNumbers.push_back(FibNumbers.at(i - 1) + FibNumbers.at(i - 2));
      }
    }
    short decoded = 0;
    for (size_t it = 0; it < chunk.size(); ++it) {
      if (chunk[it]) decoded += FibNumbers[it];
    }
    return decoded;
  }
}
