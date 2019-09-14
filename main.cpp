#include "DNAReader.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

/** The method.
* Takes in the name of a text file to read as a command-line argument. Or, if given none, prompts the user.
* After reading, continues to prompt user for text file names until told to exit.
*/
int main(int argc, char** argv) {

  DNAReader* prog;
  if (argc > 1) {
    prog = new DNAReader(argv[1]);
  } else {
    prog = new DNAReader();
  }

  while (true) {
    prog->WriteAll();
    cout << "Would you like to read through another text file? (y/n)" << endl;
    string keepGoing;
    cin >> keepGoing;
    if (tolower(keepGoing[0]) == 'y') {
      delete prog;
      DNAReader* prog = new DNAReader();
    } else {
      break;
    }
  }

}
