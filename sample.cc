#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include <map>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "Event.h"

#define FILE_EXIST 1
#define DIR_EXIST 2
#define NO_FILE_DIR_EXIST 0

using namespace std;

TFile* dataFile;
TTree* dataTree;

Event* event;
char* outDir;

vector<int> spillList;
int gl_is = 0;
int   rs;
int targ;
const char* targName[] = {"All", "LH2", "Empty", "LD2", "None", "Fe", "C", "W"};

vector<string> fileList;

bool do_analyze = false;
int fileE(const char* filename){
   struct stat st;
   if(stat(filename, &st) != 0){
      return NO_FILE_DIR_EXIST;
   }else{
      mode_t m = st.st_mode;
      if(S_ISDIR(m)){
         return NO_FILE_DIR_EXIST;
      }else{
         return FILE_EXIST;
      }
   }
}

//////////////
//   INIT   //
//////////////
void initInputRun(const char* list){
   string fileName;
   ifstream ifs;
   ifs.open(list);
   while( ifs >> fileName ) fileList.push_back(fileName);
   return;
}

void initMain(int argc, char* argv[]){
   initInputRun(argv[1]);
   outDir = argv[2];
   event = new Event();
}

bool initInputFile(int index){
   ostringstream oss;
   
   cout << " prepareing " << fileList[index] << "......" << endl;
   if( fileE(fileList[index].c_str()) == NO_FILE_DIR_EXIST ){
      cout << fileList[index] << " DOES NOT EXIST. GO TO NEXT RUN..." << endl;
      return false;
   }
   dataFile = new TFile(fileList[index].c_str(), "READ");
   dataTree = (TTree*)dataFile->Get("save");

   dataTree->SetBranchAddress("event", &event);
   cout << "  Analyze " << fileList[index].c_str() << "..." << endl;
   return true;
}

/////////////////
//   ANALYZE   //
/////////////////
void anaMain(){
   int nEvtMax = dataTree->GetEntries();
   int current =  0;
   int prev    = -1;
   for(int i = 0; i < nEvtMax ; ++i)
   {
      dataTree->GetEntry(i);
      current = (int)(double (i+1)/nEvtMax * 100);
      if( current != prev ){
         cout << "\r" << (i+1) << " / " << nEvtMax << " = " << (int)(double (i+1)/nEvtMax * 100) << " %" << flush;
         prev = (int)(double (i+1)/nEvtMax * 100);
      }

   }
   cout << endl;
}

////////////
//  MAIN  //
////////////
void deleteMemory(){
   delete   dataFile;
}

int main(int argc, char* argv[]){
   initMain(argc, argv);

   if( do_analyze ){
      for( int ir = 0 ; ir < (int)fileList.size() ; ir++ ){
         if( !initInputFile (ir) ) continue;
         anaMain();
         deleteMemory();
      }
   }   
   return EXIT_SUCCESS;
}

