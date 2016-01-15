/// std C/C++ headers
#include <stdlib.h>
#include <iostream>

/// private
#include "defaultConfigurations.h"

using namespace std;

unsigned int argumentPrintingSize = 30;

int main( int argc, char* argv[] ) {
  
  cout << endl << "[Defalut configurations:]" << endl << endl;
  for (map<string, string>::iterator it = configMap.begin();
       it!=configMap.end(); ++it){
    string argument = (*it).first;
    if (argument.size()<argumentPrintingSize){
      argument += string(argumentPrintingSize-argument.size(),' ');
    }
    cout << argument << "\t" << (*it).second << endl;
  }
  cout << endl;
  
  return 0;
}
