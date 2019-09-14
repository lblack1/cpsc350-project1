#include "DNAReader.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/** The Default Constructor.
* Prompts user for the name of a text file to read, then copies said text file into a string for easier iterating through, counting the number of lines as it goes.
*/
DNAReader::DNAReader() {

  nucleodump = "";
  ofilename = "lloydblack.out";
  totalLines = 0;
  totalA = 0;
  totalC = 0;
  totalT = 0;
  totalG = 0;
  totalAA = 0;
  totalAC = 0;
  totalAT = 0;
  totalAG = 0;
  totalCA = 0;
  totalCC = 0;
  totalCT = 0;
  totalCG = 0;
  totalTA = 0;
  totalTC = 0;
  totalTT = 0;
  totalTG = 0;
  totalGA = 0;
  totalGC = 0;
  totalGT = 0;
  totalGG = 0;
  sum = 0;
  mean = 0.0;
  variance = 0.0;

  cout << "Enter text file name: ";
  getline(cin, ifilename);

  ifstream fileIn;
  fileIn.open(ifilename);
  while (!fileIn) {
    cout << "File not found. Enter a valid text file: ";
    string ifilename;
    getline(cin, ifilename);
    fileIn.open(ifilename);
  }
  string temp;
  while (getline(fileIn, temp)) {
    nucleodump.append(temp+"\n");
    ++totalLines;
  }
  fileIn.close();

}

/** The Overloaded Constructor.
* Takes a text file as an argument, then copies said text file into a string for easier iterating through, counting the number of lines as it goes.
* param string fileName: the name of the text file to open.
*/
DNAReader::DNAReader(string fileName) {

  nucleodump = "";
  ofilename = "lloydblack.out";
  ifilename = fileName;
  totalLines = 0;
  totalA = 0;
  totalC = 0;
  totalT = 0;
  totalG = 0;
  totalAA = 0;
  totalAC = 0;
  totalAT = 0;
  totalAG = 0;
  totalCA = 0;
  totalCC = 0;
  totalCT = 0;
  totalCG = 0;
  totalTA = 0;
  totalTC = 0;
  totalTT = 0;
  totalTG = 0;
  totalGA = 0;
  totalGC = 0;
  totalGT = 0;
  totalGG = 0;
  sum = 0;
  mean = 0.0;
  variance = 0.0;

  ifstream fileIn;
  fileIn.open(ifilename);
  while (!fileIn) {
    cout << "File not found. Enter a valid text file: ";
    string ifilename;
    getline(cin, ifilename);
    fileIn.open(ifilename);
  }
  string temp;
  while (getline(fileIn, temp)) {
    nucleodump.append(temp+"\n");
    ++totalLines;
  }
  fileIn.close();

}

/** Iterates through the text file copy and counts the number of individual nucleotides and nucleotide bigrams.
* No parameters.
* Returns the sum total of all individual nucleotides in the text file. (An int)
*/
int DNAReader::ComputeTotals() {

  int tempLineCount = 0;
  for (int i = 0; i < nucleodump.length(); ++i) {

    if (tolower(nucleodump[i]) == 'a') {
      totalA++;
      if ((i+1+tempLineCount)%2 != 0 && i+1 != nucleodump.length()) { // If i is an even index and i+1 won't cause an error, essentially reads two nucleotides at a time, accounting for extra indices from line breaks
        if (tolower(nucleodump[i+1]) == 'a') {
          totalAA++;
        } else if (tolower(nucleodump[i+1]) == 'c') {
          totalAC++;
        } else if (tolower(nucleodump[i+1]) == 't') {
          totalAT++;
        } else if (tolower(nucleodump[i+1]) == 'g') {
          totalAG++;
        }
      }
    } else if (tolower(nucleodump[i]) == 'c') {
      totalC++;
      if ((i+1+tempLineCount)%2 != 0 && i+1 != nucleodump.length()) {
        if (tolower(nucleodump[i+1]) == 'a') {
          totalCA++;
        } else if (tolower(nucleodump[i+1]) == 'c') {
          totalCC++;
        } else if (tolower(nucleodump[i+1]) == 't') {
          totalCT++;
        } else if (tolower(nucleodump[i+1]) == 'g') {
          totalCG++;
        }
      }
    } else if (tolower(nucleodump[i]) == 't') {
      totalT++;
      if ((i+1+tempLineCount)%2 != 0 && i+1 != nucleodump.length()) {
        if (tolower(nucleodump[i+1]) == 'a') {
          totalTA++;
        } else if (tolower(nucleodump[i+1]) == 'c') {
          totalTC++;
        } else if (tolower(nucleodump[i+1]) == 't') {
          totalTT++;
        } else if (tolower(nucleodump[i+1]) == 'g') {
          totalTG++;
        }
      }
    } else if (tolower(nucleodump[i]) == 'g') {
      totalG++;
      if ((i+1+tempLineCount)%2 != 0 && i+1 != nucleodump.length()) {
        if (tolower(nucleodump[i+1]) == 'a') {
          totalGA++;
        } else if (tolower(nucleodump[i+1]) == 'c') {
          totalGC++;
        } else if (tolower(nucleodump[i+1]) == 't') {
          totalGT++;
        } else if (tolower(nucleodump[i+1]) == 'g') {
          totalGG++;
        }
      }
    } else if (tolower(nucleodump[i]) == '\n') {
      tempLineCount++;
    }

  }

  // Stores totals in appropriate member variables for later use.
  sum = totalA + totalC + totalT + totalG;
  totalbigrams = totalAA + totalAC + totalAT + totalAG + totalCA + totalCC + totalCT + totalCG + totalTA + totalTC + totalTT + totalTG + totalGA + totalGC + totalGT + totalGG;
  return sum;
}

/** Computes average nucleotides per line.
* No parameters.
* Returns the mean nucleotides per line. (A float)
*/
float DNAReader::ComputeMean() { // Returns the mean line length
  mean = (float)sum / (float)totalLines;
  return mean;
}

/** Computes variance of line length.
* No parameters.
* Returns the variance of nucleotides per line. (A float)
*/
float DNAReader::ComputeVariance() {

  float numerator = 0.0;
  int lineLength = 0;

  for (int i = 0; i < nucleodump.length(); ++i) {
    if (nucleodump[i] == '\n') {
      numerator += pow(((float)lineLength - mean), 2.0);
      lineLength = 0;
    } else {
      ++lineLength;
    }
  }

  variance = numerator / (float)totalLines;
  return variance;
}

/** Computes standard deviation of line length.
* No parameters.
* Returns the standard deviation of nucleotides per line. (A float)
*/
float DNAReader::ComputeStddev() {

  stddev = sqrt((double)variance);
  return stddev;

}

/** This comment applies to all ComputeRelX() and ComputeRelXY() functions.
* Computes the relative probability of encountering the corresponding nucleotide or nucleotide bigram in the text file.
* No parameters.
* Returns the relative probability of the appropriate nucleotide or nucleotide bigram. (A Float)
*/
float DNAReader::ComputeRelA() {
  relA = (float)totalA / (float)sum;
  return relA;
}

float DNAReader::ComputeRelC() {
  relC = (float)totalC / (float)sum;
  return relC;
}

float DNAReader::ComputeRelT() {
  relT = (float)totalT / (float)sum;
  return relT;
}

float DNAReader::ComputeRelG() {
  relG = (float)totalG / (float)sum;
  return relG;
}

float DNAReader::ComputeRelAA() {
  relAA = (float)totalAA / (float)totalbigrams;
  return relAA;
}

float DNAReader::ComputeRelAC() {
  relAC = (float)totalAC / (float)totalbigrams;
  return relAC;
}

float DNAReader::ComputeRelAT() {
  relAT = (float)totalAT / (float)totalbigrams;
  return relAT;
}

float DNAReader::ComputeRelAG() {
  relAG = (float)totalAG / (float)totalbigrams;
  return relAG;
}

float DNAReader::ComputeRelCA() {
  relCA = (float)totalCA / (float)totalbigrams;
  return relCA;
}

float DNAReader::ComputeRelCC() {
  relCC = (float)totalCC / (float)totalbigrams;
  return relCC;
}

float DNAReader::ComputeRelCT() {
  relCT = (float)totalCT / (float)totalbigrams;
  return relCT;
}

float DNAReader::ComputeRelCG() {
  relCG = (float)totalCG / (float)totalbigrams;
  return relCG;
}

float DNAReader::ComputeRelTA() {
  relTA = (float)totalTA / (float)totalbigrams;
  return relTA;
}

float DNAReader::ComputeRelTC() {
  relTC = (float)totalTC / (float)totalbigrams;
  return relTC;
}

float DNAReader::ComputeRelTT() {
  relTT = (float)totalTT / (float)totalbigrams;
  return relTT;
}

float DNAReader::ComputeRelTG() {
  relTG = (float)totalTG / (float)totalbigrams;
  return relTG;
}

float DNAReader::ComputeRelGA() {
  relGA = (float)totalGA / (float)totalbigrams;
  return relGA;
}

float DNAReader::ComputeRelGC() {
  relGC = (float)totalGC / (float)totalbigrams;
  return relGC;
}

float DNAReader::ComputeRelGT() {
  relGT = (float)totalGT / (float)totalbigrams;
  return relGT;
}

float DNAReader::ComputeRelGG() {
  relGG = (float)totalGG / (float)totalbigrams;
  return relGG;
}

/** Generates a random integer that fits into the Gaussian Model defined by the calculated mean and variance.
* No parameters.
* Returns a random number to be used for generating a random line of nucleotides that fits the statistical model of the text file. (an int)
*/
int DNAReader::ComputeRandomGaussian() {

  double a = (double)rand() / (double)RAND_MAX;
  double b = (double)rand() / (double)RAND_MAX;
  double c = sqrt(-2.0*log(a)) * cos(2.0*3.14159265*b);
  return sqrt(variance)*c + mean;

}

/** Generates a random line of nucleotides based on a random Gaussian line length and the relative probabilities of individual nucleotides.
* No parameters.
* Returns a string that fits the Gaussian stastical model of the text file. (a string)
*/
string DNAReader::ProduceGaussianLine() {

  string randomLine = "";
  for (int i = 0; i < ComputeRandomGaussian(); ++i) {
    double randomNucleotide = (double)rand() / (double)RAND_MAX;
    if (randomNucleotide < relA) {
      randomLine.append("A");
    } else if (randomNucleotide < relA + relC) {
      randomLine.append("C");
    } else if (randomNucleotide < relA + relC + relT) {
      randomLine.append("T");
    } else if (randomNucleotide < 1) {
      randomLine.append("G");
    }
  }
  return randomLine;

}


/** Writes all statistics and randomly generated lines to a file called "lloydblack.out".
* No parameters.
* No return.
*/
void DNAReader::WriteAll() {

  ofstream fileOut;
  fileOut.open(ofilename, ios::app);
  fileOut << " --- Lloyd Black, 2295968, CPSC350 Project 1 --- " << endl;
  fileOut << "Analysis of " << ifilename << endl; // ERROR FUCK YOU LATER LLOYD FIX THIS
  fileOut << endl;
  fileOut << "Total Nucleotides: " << this->ComputeTotals() << endl;
  fileOut << "Mean Line Length: " << this->ComputeMean() << endl;
  fileOut << "Variance: " << this->ComputeVariance() << endl;
  fileOut << "Standard Deviation: " << this->ComputeStddev() << endl;
  fileOut << endl;
  fileOut << "Relative Probability of the A Nucleotide: " << this->ComputeRelA() << endl;
  fileOut << "Relative Probability of the C Nucleotide: " << this->ComputeRelC() << endl;
  fileOut << "Relative Probability of the T Nucleotide: " << this->ComputeRelT() << endl;
  fileOut << "Relative Probability of the G Nucleotide: " << this->ComputeRelG() << endl;
  fileOut << "Relative Probability of the AA Nucleotide: " << this->ComputeRelAA() << endl;
  fileOut << "Relative Probability of the AC Nucleotide: " << this->ComputeRelAC() << endl;
  fileOut << "Relative Probability of the AT Nucleotide: " << this->ComputeRelAT() << endl;
  fileOut << "Relative Probability of the AG Nucleotide: " << this->ComputeRelAG() << endl;
  fileOut << "Relative Probability of the CA Nucleotide: " << this->ComputeRelCA() << endl;
  fileOut << "Relative Probability of the CC Nucleotide: " << this->ComputeRelCC() << endl;
  fileOut << "Relative Probability of the CT Nucleotide: " << this->ComputeRelCT() << endl;
  fileOut << "Relative Probability of the CG Nucleotide: " << this->ComputeRelCG() << endl;
  fileOut << "Relative Probability of the TA Nucleotide: " << this->ComputeRelTA() << endl;
  fileOut << "Relative Probability of the TC Nucleotide: " << this->ComputeRelTC() << endl;
  fileOut << "Relative Probability of the TT Nucleotide: " << this->ComputeRelTT() << endl;
  fileOut << "Relative Probability of the TG Nucleotide: " << this->ComputeRelTG() << endl;
  fileOut << "Relative Probability of the GA Nucleotide: " << this->ComputeRelGA() << endl;
  fileOut << "Relative Probability of the GC Nucleotide: " << this->ComputeRelGC() << endl;
  fileOut << "Relative Probability of the GT Nucleotide: " << this->ComputeRelGT() << endl;
  fileOut << "Relative Probability of the GG Nucleotide: " << this->ComputeRelGG() << endl;
  fileOut << endl;
  fileOut << " --- Normally Modeled Strings --- " << endl;
  for (int i = 0; i < 1000; ++i) {
    fileOut << this->ProduceGaussianLine() << endl;
  }
  fileOut << endl;
  fileOut << endl;
  fileOut << endl;

  fileOut.close();

}
