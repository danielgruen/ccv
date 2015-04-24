#include "filter.h"
#include <CCfits/CCfits>

using namespace std;
using namespace CCfits;

int max(int a, int b)
{
 if(a>=b) return a;
 return b;
}


double pow2_filter(double x)
{
 return x*x; 
}

// ObjectPrototype

  ObjectPrototype::ObjectPrototype(const int keys[], const string names[], const Type types[], int n, string prefix) : intSortVkey(-1), doubleSortVkey(-1) 
  {
   for (int i=0; i<n; i++)
   {
    switch (types[i]) 
    {
      case integerType:
	cerr << "# reading integer column " << keys[i] << " as " << prefix+names[i] << endl;
	intPropertyName.push_back(prefix+names[i]);
	intPropertyKey.push_back(keys[i]-1);
	break;
      case doubleType:
	cerr << "# reading double column " << keys[i] << " as " << prefix+names[i] << endl;
	doublePropertyName.push_back(prefix+names[i]);
	doublePropertyKey.push_back(keys[i]-1);
	break;
      case stringType:
	cerr << "# reading string column " << keys[i] << " as " << prefix+names[i] << endl;
	stringPropertyName.push_back(prefix+names[i]);
	stringPropertyKey.push_back(keys[i]-1);
	break;
    }
   } 
  }

  int ObjectPrototype::intVkey(string name, bool append)
  {
   for (int i=0; i<intPropertyName.size(); i++)
   {
    if (intPropertyName[i]==name) return i;
   }
   if(append)
    {
     intPropertyName.push_back(name);
     intPropertyKey.push_back(-1);
    }
    return -1; // not found
  }  
  int ObjectPrototype::doubleVkey(string name, bool append) 
  {
   //cerr << "checking double key " << name << endl;
   for (int i=0; i<doublePropertyName.size(); i++)
   {
    //cerr << "is it like " << doublePropertyName[i] << "?" << endl;
    if (doublePropertyName[i]==name) return i;
    //cerr << "no!" << endl;
   }
   if(append)
    {
     doublePropertyName.push_back(name);
     doublePropertyKey.push_back(-1);
    }
    //cerr << "not found" << endl;
    return -1; // not found
  }  
  int ObjectPrototype::stringVkey(string name, bool append)
  {
   for (int i=0; i<stringPropertyName.size(); i++)
   {
    if (stringPropertyName[i]==name) return i;
   }
   if(append)
    {
     stringPropertyName.push_back(name);
     stringPropertyKey.push_back(-1);
    }
    return -1; // not found
  }

// Object

  Object::Object(ObjectPrototype *p) : prototype(p) 
  { // create empty object according to a prototype
      intProperty.resize(prototype->intPropertyKey.size(),0);
      doubleProperty.resize(prototype->doublePropertyKey.size(),0.);
      stringProperty.resize(prototype->stringPropertyKey.size(),"");
  }  


  Object::Object(ObjectPrototype *p, ifstream &file)
  {
   prototype=p;
   readLine(file);
  }
  
  int Object::intPropertyValue(int vkey) const
  {
   if(vkey>=0 && vkey<intProperty.size())
   return intProperty[vkey];
   return 0./0.;
  } 
  int Object::intPropertyValue(string name) const
  {
   return intPropertyValue(prototype->intVkey(name)); 
  }
    
  double Object::doublePropertyValue(int vkey) const
  {
   if(vkey>=0 && vkey<doubleProperty.size())
   return doubleProperty[vkey];
   return 0./0.;
  } 
  double Object::doublePropertyValue(string name) const
  {
   return doublePropertyValue(prototype->doubleVkey(name)); 
  }
    
  string Object::stringPropertyValue(int vkey) const
  {
   if(vkey>=0 && vkey<stringProperty.size())
   return stringProperty[vkey]; 
   return "";
  } 
  string Object::stringPropertyValue(string name) const
  {
   return stringPropertyValue(prototype->stringVkey(name)); 
  }
  double Object::numericValue(string name) const // convenience wrapper
  {
    int ik=prototype->intVkey(name);
    if(ik>=0) return intPropertyValue(ik);
    int dk=prototype->doubleVkey(name);
    if(dk>=0) return doublePropertyValue(dk);
    return 0./0.;
  }
  
  
  
  bool Object::setIntProperty(int vkey, int value){
   if(vkey>=0 && vkey<intProperty.size())
   {
    intProperty[vkey]=value; return true; 
   }
   return false;
  }
  bool Object::setIntProperty(string name, int value){
   return setIntProperty(prototype->intVkey(name), value); 
  }  
  
  bool Object::setDoubleProperty(int vkey, double value){
   if(vkey>=0 && vkey<doubleProperty.size())
   {
    doubleProperty[vkey]=value; return true; 
   }
   return false;
  }
  bool Object::setDoubleProperty(string name, double value){
   return setDoubleProperty(prototype->doubleVkey(name), value); 
  }  
  
  bool Object::setStringProperty(int vkey, string value){
   if(vkey>=0 && vkey<stringProperty.size())
   {
    stringProperty[vkey]=value; return true; 
   }
   return false;
  }
  bool Object::setStringProperty(string name, string value){
   return setStringProperty(prototype->stringVkey(name), value); 
  }
  
  
  void Object::readLine(ifstream &file)
  {
    string line;
    getlineNoComment(file, line); 
    readObject(line);
  }
  
  void Object::makeEmpty()
  {
        intProperty.clear(); doubleProperty.clear(); stringProperty.clear();
    
     	for(vector<int>::iterator it = prototype->intPropertyKey.begin(); it != prototype->intPropertyKey.end(); ++it) {
	  intProperty.push_back(0);
        }
     	for(vector<int>::iterator it = prototype->doublePropertyKey.begin(); it != prototype->doublePropertyKey.end(); ++it) {
	  doubleProperty.push_back(0.);
        }
     	for(vector<int>::iterator it = prototype->stringPropertyKey.begin(); it != prototype->stringPropertyKey.end(); ++it) {
	  stringProperty.push_back("");
        }
  }
  
  void Object::readObject(string line)
  // read full line with stringstuff::getlineNoComment(ifstream file, line), then call this
  {
      istringstream iss(line);
      int nread=0;
      for(std::vector<int>::iterator it = prototype->intPropertyKey.begin(); it != prototype->intPropertyKey.end(); ++it) {
	nread=max(nread,*it);
      }
      for(std::vector<int>::iterator it = prototype->doublePropertyKey.begin(); it != prototype->doublePropertyKey.end(); ++it) {
	nread=max(nread,*it);
      }
      for(std::vector<int>::iterator it = prototype->stringPropertyKey.begin(); it != prototype->stringPropertyKey.end(); ++it) {
	nread=max(nread,*it);
      }

      intProperty.resize(prototype->intPropertyKey.size());
      doubleProperty.resize(prototype->doublePropertyKey.size());
      stringProperty.resize(prototype->stringPropertyKey.size());

      for (int i=0; i<=nread; i++)
      {
	string buf;
	iss >> buf;
	
	int key=0;
	for(vector<int>::iterator it = prototype->intPropertyKey.begin(); it != prototype->intPropertyKey.end(); ++it) {
	  if(*it==i) { intProperty[key]=atoi(buf.c_str()); }
          key++;
        }
        key=0;
	for(vector<int>::iterator it = prototype->doublePropertyKey.begin(); it != prototype->doublePropertyKey.end(); ++it) {
	  if(*it==i) { doubleProperty[key]=atof(buf.c_str()); }
          key++;
        }
        key=0;
	for(vector<int>::iterator it = prototype->stringPropertyKey.begin(); it != prototype->stringPropertyKey.end(); ++it) {
	  if(*it==i) { stringProperty[key]=buf; }
          key++;
        }
      }
  }  

  ObjectCollection::ObjectCollection(ifstream &file, const int keys[], const string names[], const Type types[], int n, string prefix)
  {
   if(!file) {
    cerr << "file not open! exiting" << endl; exit(1); 
   }
   prototype= new ObjectPrototype(keys,names,types,n,prefix);
   my_memory=true;
   appendFile(file);
  }
  
  ObjectCollection::ObjectCollection(string file, string extension, vector<string> columns, string prefix)
  {
   vector<string> hdus(1);
   hdus[0]=extension;
   auto_ptr<FITS> pInfile(new FITS(file, Read, hdus, false));
   ExtHDU& table = pInfile->extension(extension);
   
   
   
   map<string, Column *> colmap = table.column();
   
   prototype=new ObjectPrototype();
   
   vector<vector<int> > intcolumns;
   vector<vector<double> > doublecolumns;
   vector<vector<string> > stringcolumns;
   
   // build prototype
   int extracolumns=0;
   for (map<string, Column *>::iterator iter = colmap.begin(); 
	(iter != colmap.end() && (columns.size()==0 || intcolumns.size()+doublecolumns.size()+stringcolumns.size()-extracolumns<columns.size()));
        iter++) {

    //cerr << "# found column " << (*iter).first << endl;
     
    if (columns.size()>0 && std::find(columns.begin(), columns.end(), (*iter).first) == columns.end()) // do not care about this column
      continue;
    
    cerr << "# reading column " << (*iter).first << endl;
    cerr << "# type " << (*iter).second->type() << endl;
    switch((*iter).second->type()) {
      case Tbit:
      case Tbyte:
      case Tlogical:
      case Tushort:
      case Tshort:
      case Tuint:
      case Tint:
      case Tulong:
      case Tlong:
      case Tlonglong:
      {
	prototype->intVkey(prefix+(*iter).first,true);
	vector<int> v;
	intcolumns.push_back(v);
	
	try {
	(*iter).second->read(intcolumns.back(),1,table.rows());
	} catch (...) {
	  cout << "# caught an error reading int column, maybe it's an array column" << endl; 
	  vector<valarray<int> > vals; // vector of array cells
	  (*iter).second->readArrays(vals,1,table.rows());
	  vector<int> v1[vals[0].size()-1];
	  for (int i=0; i<vals.size(); i++)
	  {
	    (intcolumns.back()).push_back(vals[i][0]);
	    for(int j=0; j<vals[0].size()-1; j++)
	    {
	    v1[j].push_back(vals[i][j+1]); 
	    }
	  }
	  string number[] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" };
	  if(vals[0].size()>15) {
	   cout << "too many columns" << endl; exit(1); 
	  }
	  for(int j=0; j<vals[0].size()-1; j++)
	  {
		prototype->intVkey(prefix+(*iter).first+number[j],true);
		intcolumns.push_back(v1[j]);
		cerr << "# reading column " << (*iter).first+number[j] << endl;
		cerr << "# size: " << v1[j].size() << endl;
		extracolumns++;
	  }
	}
      }
	break;	
      case Tfloat:
      case Tdouble:
      {
	prototype->doubleVkey(prefix+(*iter).first,true);
	vector<double> v;
	doublecolumns.push_back(v);
	try {
	(*iter).second->read(doublecolumns.back(),1,table.rows());
	} catch (...) {
	  cout << "# caught an error reading double column, maybe it's an array column" << endl; 
	  vector<valarray<double> > vals; // vector of array cells
	  (*iter).second->readArrays(vals,1,table.rows());
	  vector<double> v1[vals[0].size()-1];
	  for (int i=0; i<vals.size(); i++)
	  {
	    (doublecolumns.back()).push_back(vals[i][0]);
	    for(int j=0; j<vals[0].size()-1; j++)
	    {
	    v1[j].push_back(vals[i][j+1]); 
	    }
	  }
	  string number[] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" };
	  if(vals[0].size()>15) {
	   cout << "too many columns" << endl; exit(1); 
	  }
	  for(int j=0; j<vals[0].size()-1; j++)
	  {
		prototype->doubleVkey(prefix+(*iter).first+number[j],true);
		doublecolumns.push_back(v1[j]);
		cerr << "# reading column " << (*iter).first+number[j] << endl;
		extracolumns++;
	  }
	}
      }
	break;	
      case Tstring:
      {
	prototype->stringVkey(prefix+(*iter).first,true);
	vector<string> v;
	stringcolumns.push_back(v);
	(*iter).second->read(stringcolumns.back(),1,table.rows());
      }
	break;	
      default:
	cerr << "cannot handle column type " << (*iter).second->type() << "; exiting..." << endl; exit(1);
    }
    
    //cerr << *prototype << endl;
   }
   
   readFITSTable(table.rows(), intcolumns, doublecolumns, stringcolumns);
   my_memory=true;
  }
    
  ObjectCollection::ObjectCollection(string file, int extension, vector<string> columns, string prefix)
  {
   auto_ptr<FITS> pInfile(new FITS(file, Read, extension, false));
   ExtHDU& table = pInfile->extension(extension);
   
   map<string, Column *> colmap = table.column();
   
   prototype=new ObjectPrototype();
   
   vector<vector<int> > intcolumns;
   vector<vector<double> > doublecolumns;
   vector<vector<string> > stringcolumns;
   
   // build prototype
   int extracolumns=0;
   for (map<string, Column *>::iterator iter = colmap.begin(); 
	(iter != colmap.end() && (columns.size()==0 || intcolumns.size()+doublecolumns.size()+stringcolumns.size()-extracolumns<columns.size()));
        iter++) {

    //cerr << "# found column " << (*iter).first << endl;
     
    if (columns.size()>0 && std::find(columns.begin(), columns.end(), (*iter).first) == columns.end()) // do not care about this column
      continue;
    
    cerr << "# reading column " << (*iter).first << endl;
    //cerr << "# type " << (*iter).second->type() << endl;
    switch((*iter).second->type()) {
      case Tbit:
      case Tbyte:
      case Tlogical:
      case Tushort:
      case Tshort:
      case Tuint:
      case Tint:
      case Tulong:
      case Tlong:
      case Tlonglong:
      {
	prototype->intVkey(prefix+(*iter).first,true);
	vector<int> v;
	intcolumns.push_back(v);
	
	try {
	(*iter).second->read(intcolumns.back(),1,table.rows());
	} catch (...) {
	  cout << "# caught an error reading int column, maybe it's an array column" << endl; 
	  vector<valarray<int> > vals; // vector of array cells
	  (*iter).second->readArrays(vals,1,table.rows());
	  vector<int> v1[vals[0].size()-1];
	  for (int i=0; i<vals.size(); i++)
	  {
	    (intcolumns.back()).push_back(vals[i][0]);
	    for(int j=0; j<vals[0].size()-1; j++)
	    {
	    v1[j].push_back(vals[i][j+1]); 
	    }
	  }
	  string number[] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" };
	  if(vals[0].size()>15) {
	   cout << "too many columns" << endl; exit(1); 
	  }
	  for(int j=0; j<vals[0].size()-1; j++)
	  {
		prototype->intVkey(prefix+(*iter).first+number[j],true);
		intcolumns.push_back(v1[j]);
		cerr << "# reading column " << (*iter).first+number[j] << endl;
		cerr << "# size: " << v1[j].size() << endl;
		extracolumns++;
	  }
	}
      }
	break;	
      case Tfloat:
      case Tdouble:
      {
	prototype->doubleVkey(prefix+(*iter).first,true);
	vector<double> v;
	doublecolumns.push_back(v);
	try {
	(*iter).second->read(doublecolumns.back(),1,table.rows());
	} catch (...) {
	  cout << "# caught an error reading double column, maybe it's an array column" << endl; 
	  vector<valarray<double> > vals; // vector of array cells
	  (*iter).second->readArrays(vals,1,table.rows());
	  vector<double> v1[vals[0].size()-1];
	  for (int i=0; i<vals.size(); i++)
	  {
	    (doublecolumns.back()).push_back(vals[i][0]);
	    for(int j=0; j<vals[0].size()-1; j++)
	    {
	    v1[j].push_back(vals[i][j+1]); 
	    }
	  }
	  string number[] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" };
	  if(vals[0].size()>15) {
	   cout << "too many columns" << endl; exit(1); 
	  }
	  for(int j=0; j<vals[0].size()-1; j++)
	  {
		prototype->doubleVkey(prefix+(*iter).first+number[j],true);
		doublecolumns.push_back(v1[j]);
		cerr << "# reading column " << (*iter).first+number[j] << endl;
		extracolumns++;
	  }
	}
      }
	break;	
      case Tstring:
      {
	prototype->stringVkey(prefix+(*iter).first,true);
	vector<string> v;
	stringcolumns.push_back(v);
	(*iter).second->read(stringcolumns.back(),1,table.rows());
      }
	break;	
      default:
	cerr << "cannot handle column type " << (*iter).second->type() << "; exiting..." << endl; exit(1);
    }
    
    //cerr << *prototype << endl;
   }
   
   readFITSTable(table.rows(), intcolumns, doublecolumns, stringcolumns);
   my_memory=true;
  }
  
  void ObjectCollection::readFITSTable(int rows, vector<vector<int> > &intcolumns,  vector<vector<double> > &doublecolumns,  vector<vector<string> > &stringcolumns)
  {  
   // build objects and read data
   for (int i=1; i<=rows; i++)
   {
     Object *o=new Object(prototype);
     //o->makeEmpty();
     
     for(int it=0; it<intcolumns.size(); it++){ 
       //cout << "# int column " << it << ": size " << intcolumns[it].size() << endl;
       o->setIntProperty(it,(intcolumns[it])[i-1]);
     }
          
     for(int it=0; it<doublecolumns.size(); it++){ 
       //cout << "# double column " << it << ": size " << doublecolumns[it].size() << endl;
       o->setDoubleProperty(it,(doublecolumns[it])[i-1]);
     }
          
     for(int it=0; it<stringcolumns.size(); it++){ 
       o->setStringProperty(it,(stringcolumns[it])[i-1]);
     }

     objects.push_back(o);
   } 
  }
  
  void ObjectCollection::appendFITSTable(string file, string extension, vector<string> columns)
  {
   vector<string> hdus(1);
   hdus[0]=extension;
   auto_ptr<FITS> pInfile(new FITS(file, Read, hdus, false));
   ExtHDU& table = pInfile->extension(extension);
   map<string, Column *> colmap = table.column();
   
   vector<vector<int> > intcolumns;
   vector<vector<double> > doublecolumns;
   vector<vector<string> > stringcolumns;
   
   for (map<string, Column *>::iterator iter = colmap.begin(); (iter != colmap.end() && (columns.size()==0 || intcolumns.size()+doublecolumns.size()+stringcolumns.size()<columns.size())); iter++) {

    if (columns.size()>0 && std::find(columns.begin(), columns.end(), (*iter).first) == columns.end()) // do not care about this column
      continue;
    
    cerr << "# reading column " << (*iter).first << endl;
    
    switch((*iter).second->type()) {
      case Tbit:
      case Tbyte:
      case Tlogical:
      case Tushort:
      case Tshort:
      case Tuint:
      case Tint:
      case Tulong:
      case Tlong:
      case Tlonglong:
      {
	vector<int> v;
	intcolumns.push_back(v);
	(*iter).second->read(intcolumns.back(),1,table.rows());
      }
	break;	
      case Tfloat:
      case Tdouble:
      {
	vector<double> v;
	doublecolumns.push_back(v);
	(*iter).second->read(doublecolumns.back(),1,table.rows());
      }
	break;	
      case Tstring:
      {
	vector<string> v;
	stringcolumns.push_back(v);
	(*iter).second->read(stringcolumns.back(),1,table.rows());
      }
	break;	
      default:
	cerr << "cannot handle column type " << (*iter).second->type() << "; exiting..." << endl; exit(1);
    }
   }
   
   readFITSTable(table.rows(), intcolumns, doublecolumns, stringcolumns);
  }
  
  ObjectCollection::~ObjectCollection()
  {
   if(!my_memory) return; // only delete if prototype+objects are not inherited
   delete prototype;
   for(vector<Object*>::iterator it = objects.begin(); it != objects.end(); ++it) {
	delete *it; 
   }
  }
  
  void ObjectCollection::setPrototype(const int keys[], const string names[], const Type types[], int n)
  {
   if(!prototype) prototype=new ObjectPrototype(keys, names, types, n);  
  }
  
  void ObjectCollection::appendFile(ifstream &file)
  {
   if(!file) {
    cerr << "file not open! exiting" << endl; exit(1); 
   }
   while(1)
   {
    Object *o = new Object(prototype, file);
    if(!file.good()) return;
    objects.push_back(o);
   }
  }
  
  void ObjectCollection::appendPrototype(ObjectPrototype *otherPrototype, string prefix, string exception) {
    // appends prototype and creates columns for new fields
    
    	for(vector<string>::iterator it = otherPrototype->intPropertyName.begin(); it != otherPrototype->intPropertyName.end(); ++it) {
	  if((*it)!=exception)
	    createIntPropertyIfNecessary(prefix+(*it));
        }
        for(vector<string>::iterator it = otherPrototype->doublePropertyName.begin(); it != otherPrototype->doublePropertyName.end(); ++it) {
	  if((*it)!=exception)
	    createDoublePropertyIfNecessary(prefix+(*it));
        }
        for(vector<string>::iterator it = otherPrototype->stringPropertyName.begin(); it != otherPrototype->stringPropertyName.end(); ++it) {
	  if((*it)!=exception)
	    createStringPropertyIfNecessary(prefix+(*it));
        }
  }
  
  void ObjectCollection::extendCollection(ObjectCollection &other, string prefix) {
    
	if(other.size() != this->size()) {
	 cerr << "cannot extend catalog by columns of different row count; ignoring" << endl;
	 return;
	}

        for(vector<string>::iterator it = other.prototype->intPropertyName.begin(); it != other.prototype->intPropertyName.end(); ++it) {
	  int c = createIntPropertyIfNecessary(prefix+(*it));
	  int cother = other.prototype->intVkey(*it);
	  for(int i=0; i<objects.size(); i++) 
	  {
	    objects[i]->setIntProperty(c,other.objects[i]->intPropertyValue(cother));
	  }
        }
        for(vector<string>::iterator it = other.prototype->doublePropertyName.begin(); it != other.prototype->doublePropertyName.end(); ++it) {
	  int c = createDoublePropertyIfNecessary(prefix+(*it));
	  int cother = other.prototype->doubleVkey(*it);
	  for(int i=0; i<objects.size(); i++) 
	  {
	    objects[i]->setDoubleProperty(c,other.objects[i]->doublePropertyValue(cother));
	  }
        }
        for(vector<string>::iterator it = other.prototype->stringPropertyName.begin(); it != other.prototype->stringPropertyName.end(); ++it) {
	  int c = createStringPropertyIfNecessary(prefix+(*it));
	  int cother = other.prototype->stringVkey(*it);
	  for(int i=0; i<objects.size(); i++) 
	  {
	    objects[i]->setStringProperty(c,other.objects[i]->stringPropertyValue(cother));
	  } 
        }
	
  }
 
  void ObjectCollection::dualMatchCollection(ObjectCollection &other, string prefix, string xname, string yname, string xnamenew, string ynamenew, double r1, double r2, bool append_unmatched, bool shout)
  {
  
       // add information from additional catalog to columns with prefix "prefix"
       //  match to current objects if position xname,yname agrees within less than r1 and within less than r2 uniquely
       //  otherwise append catalog
       
       // check whether match key exists
       int selfmatch=createDoublePropertyIfNecessary("match", 1, true);
       
       // generate key indicating there is a match in the added catalog
       int othermatch=createDoublePropertyIfNecessary(prefix+"match", 0, true);
       
        
       //if(shout)
       cerr << "# other catalog " << prefix << " has " << other.size() << " entries, I have " << this->size() << endl;
       
       appendPrototype(other.prototype, prefix, "match"); // do not add an extra match column if the matched catalog has one...
  
       int myx=prototype->doubleVkey(xname);
       int myy=prototype->doubleVkey(yname);
       int otherx=other.prototype->doubleVkey(xnamenew);
       int othery=other.prototype->doubleVkey(ynamenew);
       
       //if(shout)
       cerr << "# matching double keys " << myx << "," << myy << " from me with " << otherx << "," << othery << " from other catalog" << endl;
       
       double rmax = r2;

       r1=r1*r1*1.5;
       r2=r2*r2*1.5;
       
       int nmatch=0;
       int nnomatch=0;

       sort(xname);
       
       //for(vector<Object*>::iterator it = other.objects.begin(); it != other.objects.end(); ++it) // actually match/append other catalog entries
       
       int mul=objects.size()-1;
       
#pragma omp parallel for
       for(int iit=0; iit<other.size(); iit++)
       { 
/*	 if(shout) {
#pragma omp critical
{
	  for(int iii=0; iii<size(); iii++)
	  {
	   if(objects[iii]->doublePropertyValue(othermatch)<0) {
	    cerr << "found object " << iii/size() << " at othermatch<0" << endl;
	    cerr << *(objects[iii]) << endl;
	    exit(0);
	   }
	  }
}
	 }*/
	 
	 Object *it=other[iit];
	 Object *match=0;
         // find id range for iteration
         int llid=0; // lower limit of starting point
         int ulid=mul; // upper limit of starting point
         int uuid=mul; // upper limit of end point
         double ox = (it)->doublePropertyValue(otherx);
         double oy = (it)->doublePropertyValue(othery);
         while(ulid-llid>2) {
	   int tlid=(llid+ulid)/2;
           double tx=objects[tlid]->doublePropertyValue(myx)-ox+rmax;
           if(tx<=0) llid=tlid;
           else ulid=tlid;
           //cerr << "# lid between " << llid << " and " << ulid << endl;
         }
	 int luid = llid+1;
         while(uuid-luid>2) {
	   int tuid=(luid+uuid)/2;
           double tx=objects[tuid]->doublePropertyValue(myx)-ox-rmax;
           if(tx<=0) luid=tuid;
           else uuid=tuid;
           //cerr << "# uid between " << luid << " and " << uuid << endl;
         }


	 for(int imyit = llid; imyit <= uuid; ++imyit)
	 {
	  double R2 = pow2_filter(ox-(objects[imyit])->doublePropertyValue(myx))+pow2_filter(oy-(objects[imyit])->doublePropertyValue(myy));

	  if(shout)
	  cerr << "# considering object with matchpointer=" << match << " R2=" << R2 << " selfmatch=" <<  objects[imyit]->doublePropertyValue(selfmatch) << " othermatch=" << objects[imyit]->doublePropertyValue(othermatch) << endl;
	  
	  if(R2<r1 && !match && !objects[imyit]->doublePropertyValue(othermatch)) { // new match!; added !objects[imyit]->doublePropertyValue(othermatch) on Nov 04 2013
	    match=(objects[imyit]);
	  } else if (R2<r2) { // second match or bad match!
	    if(shout) {
	      cerr << "# found second or bad match, first match is " << match << ", radius is " << R2 << ", (X,Y) is " << objects[imyit]->doublePropertyValue(myx) << "," << objects[imyit]->doublePropertyValue(myy) << " in me and " << (it)->doublePropertyValue(otherx) << "," << (it)->doublePropertyValue(othery) << " in other catalog; othermatch is " << objects[imyit]->doublePropertyValue(othermatch) << endl;
	    }
	    match=0;
	    break;
	  }
	 }
	 
	 if(match) {
#pragma omp critical
	  {
	  ++nmatch;
	  if(shout)
	    cerr << "# found match in me for object at " << (it)->doublePropertyValue(otherx) << "," << (it)->doublePropertyValue(othery) << endl;
	  }
	  for(vector<string>::iterator otherit = other.prototype->intPropertyName.begin(); otherit != other.prototype->intPropertyName.end(); ++otherit) {
	    match->setIntProperty(prefix+*otherit,(it)->intPropertyValue(*otherit));
          }
          for(vector<string>::iterator otherit = other.prototype->doublePropertyName.begin(); otherit != other.prototype->doublePropertyName.end(); ++otherit) {
	    match->setDoubleProperty(prefix+*otherit,(it)->doublePropertyValue(*otherit));
          }
          for(vector<string>::iterator otherit = other.prototype->stringPropertyName.begin(); otherit != other.prototype->stringPropertyName.end(); ++otherit) {
	    match->setStringProperty(prefix+*otherit,(it)->stringPropertyValue(*otherit));
          }
          match->setDoubleProperty(othermatch,1);
	 }
	 else {
#pragma omp critical
	   {
	   nnomatch++;
	   if(shout)
	     cerr << "# found no match in me for object at " << (it)->doublePropertyValue(otherx) << "," << (it)->doublePropertyValue(othery) << endl;
	   }
	   if(append_unmatched) {
	   Object *o;
#pragma omp critical
	   o = appendEmpty();
	   for(vector<string>::iterator otherit = other.prototype->intPropertyName.begin(); otherit != other.prototype->intPropertyName.end(); ++otherit) {
	    o->setIntProperty(prefix+*otherit,(it)->intPropertyValue(*otherit));
           }
           for(vector<string>::iterator otherit = other.prototype->doublePropertyName.begin(); otherit != other.prototype->doublePropertyName.end(); ++otherit) {
	    o->setDoubleProperty(prefix+*otherit,(it)->doublePropertyValue(*otherit));
           }
           for(vector<string>::iterator otherit = other.prototype->stringPropertyName.begin(); otherit != other.prototype->stringPropertyName.end(); ++otherit) {
	    o->setStringProperty(prefix+*otherit,(it)->stringPropertyValue(*otherit));
           }
	   o->setDoubleProperty(othermatch,1);
	   o->setDoubleProperty(selfmatch,0);
	   }
	 }
	 
       }
       
       cerr << "# found " << nmatch << " matching and " << nnomatch << " unmatching objects" << endl;
  

  } 
 
  void ObjectCollection::dualMatchFile(ifstream &file, string prefix, const int keys[], const string names[], const Type types[], int n, string xname, string yname,  string xnamenew, string ynamenew, double r1, double r2, bool append_unmatched, bool shout)
  {
       if(!file) {
	cerr << "file not open! exiting" << endl; exit(1); 
       }
      
       ObjectCollection other(file, keys, names, types, n);

       dualMatchCollection(other, prefix, xname, yname, xnamenew, ynamenew, r1, r2, append_unmatched, shout);
  }
    
  void ObjectCollection::appendCollectionDeep(ObjectCollection &cc)
  {
   for(int i=0; i<cc.size(); i++) {
    Object *o = new Object(*(cc[i]));
    o->prototype = prototype; // tacitly assumes the have structurally the same prototype, even though it might be two separate copies of it
    objects.push_back(o);
   }
  }
  
  void ObjectCollection::appendObject(Object *cc)
  {
   objects.push_back(cc); 
  }
  
  void ObjectCollection::makeEmpty(int nrows)
  {
   if(my_memory)
   {
     for(vector<Object*>::iterator it = objects.begin(); it != objects.end(); ++it)
	delete *it; 
   }
   objects.clear();
   for(int i=0; i<nrows; i++)
   {
    appendEmpty();
   }
   if(nrows>0) my_memory=true;
  }
  
  Object* ObjectCollection::appendEmpty() 
  {
    Object *o = new Object(prototype);
    //o->makeEmpty();
    objects.push_back(o);
    return o;
  }
  
  /////////////////////////////////// FILTERING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
  ObjectCollection* ObjectCollection::filter(Filter f) const
  {
   int key = prototype->intVkey(f.name);
   if(key<0) 
   {
     key=prototype->doubleVkey(f.name);
     return filterDouble(key, f.minVal, f.maxVal); 
   }
   return filterInt(key, int(f.minVal+0.999), int(f.maxVal));   
  }
  
  ObjectCollection* ObjectCollection::filter(vector<Filter> f) const
  {
   vector<Filter> doubleF;
   vector<int> doubleVkeys;
   vector<Filter> intF;
   vector<int> intVkeys;
   for(int i=0; i<f.size(); i++)
   {
    int key=prototype->doubleVkey(f[i].name);
    if(key>=0) 
    {
      doubleF.push_back(f[i]); // it's a double filter
      doubleVkeys.push_back(key);
      continue;
    }
    key=prototype->intVkey(f[i].name);
    if(key>=0) 
    {
      intF.push_back(f[i]); // it's a double filter
      intVkeys.push_back(key);
      continue;
    }
    
    cerr << "filtering by unknown column " << f[i].name << ", aborting!" << endl; return 0;
   }
   
   // now filter
   ObjectCollection* o = new ObjectCollection();
   o->prototype = prototype;
   for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
     for(int fit=0; fit<doubleF.size(); fit++)
     {
       if((*it)->doubleProperty[doubleVkeys[fit]] < doubleF[fit].minVal || (*it)->doubleProperty[doubleVkeys[fit]] > doubleF[fit].maxVal) // deselect!
	goto fd;
     }
     for(int fit=0; fit<intF.size(); fit++)
     {
       if((*it)->intProperty[intVkeys[fit]] < intF[fit].minVal || (*it)->intProperty[intVkeys[fit]] > intF[fit].maxVal) // deselect!
	goto fd;
     }     
     
     o->objects.push_back(*it); // if not deselected, keep

     fd: 
     ;
   }
   return o;
  }
  
    
  vector<unsigned char> ObjectCollection::filterFlagVector(vector<Filter> f) const // returns vector which is 1 for good objects, 0 for deselected objects
  {
   vector<unsigned char> v;
    
   vector<Filter> doubleF;
   vector<int> doubleVkeys;
   vector<Filter> intF;
   vector<int> intVkeys;
   for(int i=0; i<f.size(); i++)
   {
    int key=prototype->doubleVkey(f[i].name);
    if(key>=0) 
    {
      doubleF.push_back(f[i]); // it's a double filter
      doubleVkeys.push_back(key);
      continue;
    }
    key=prototype->intVkey(f[i].name);
    if(key>=0) 
    {
      intF.push_back(f[i]); // it's a double filter
      intVkeys.push_back(key);
      continue;
    }
   }
   
   for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
     for(int fit=0; fit<doubleF.size(); fit++)
     {
       if((*it)->doubleProperty[doubleVkeys[fit]] < doubleF[fit].minVal || (*it)->doubleProperty[doubleVkeys[fit]] > doubleF[fit].maxVal) // deselect!
       {
	v.push_back(false);
	goto fd;
       }
     }
     for(int fit=0; fit<intF.size(); fit++)
     {
       if((*it)->intProperty[intVkeys[fit]] < intF[fit].minVal || (*it)->intProperty[intVkeys[fit]] > intF[fit].maxVal) // deselect!
       {
	v.push_back(false);
	goto fd;
       }
     }     
     
     v.push_back(true); // if not deselected, it's good

     fd: 
     ;
   }
   return v;
  }
  
  ObjectCollection* ObjectCollection::filterInt(string filterkey, int minInt, int maxInt) const
  {
   return filterInt(prototype->intVkey(filterkey), minInt, maxInt); 
  }
  ObjectCollection* ObjectCollection::filterInt(int filterkey, int minInt, int maxInt) const
  {
   if(filterkey<0) {
        cerr << "filtering for invalid int key between " << minInt << " and " << maxInt << endl;
	return 0; // invalid name
   }
   ObjectCollection* o = new ObjectCollection();
   o->prototype = prototype; // same prototype
   for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	if ((*it)->intProperty[filterkey]>=minInt && (*it)->intProperty[filterkey]<=maxInt)
	  o->objects.push_back(*it);
   }
   return o; 
  }
  
  ObjectCollection* ObjectCollection::filterDouble(string filterkey, double minDouble, double maxDouble) const
  {
   return filterDouble(prototype->doubleVkey(filterkey), minDouble, maxDouble); 
  }
  ObjectCollection* ObjectCollection::filterDouble(int filterkey, double minDouble, double maxDouble) const
  {
   if(filterkey<0) {
	cerr << "filtering for invalid double key between " << minDouble << " and " << maxDouble << endl;
	return 0; // invalid name
   }
   ObjectCollection* o = new ObjectCollection();
   o->prototype = prototype; // same prototype
   for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	if ((*it)->doubleProperty[filterkey]>=minDouble && (*it)->doubleProperty[filterkey]<=maxDouble)
	  o->objects.push_back(*it);
   }
   return o; 
  }

  ObjectCollection* ObjectCollection::filterString(string filterkey, string filtervalue) const
  {
   return filterString(prototype->stringVkey(filterkey), filtervalue);
  }
  ObjectCollection* ObjectCollection::filterString(int filterkey, string filtervalue) const
  {
   if(filterkey<0) {
	cerr << "filtering for invalid string key" << endl;
        return 0;
   }
   ObjectCollection* o = new ObjectCollection();
   o->prototype = prototype; // same prototype
   for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	if ((*it)->stringProperty[filterkey]==filtervalue)
	  o->objects.push_back(*it);
   }
   return o; 
  }  
  
  ObjectCollection* ObjectCollection::filterStringRemove(string filterkey, vector<string> baddata) const
  {
   return filterStringRemove(prototype->stringVkey(filterkey), baddata);
  }
  ObjectCollection* ObjectCollection::filterStringRemove(int filterkey, vector<string> baddata) const
  {
   if(filterkey<0) {
	cerr << "filtering for invalid string key" << endl;
        return 0;
   }
   ObjectCollection* o = new ObjectCollection();
   o->prototype = prototype; // same prototype
   
   for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
     string data = (*it)->stringProperty[filterkey];
     bool found=false;
     for(vector<string>::const_iterator is = baddata.begin(); is != baddata.end(); ++is) {
      if (found || data.find((*is)) != std::string::npos) {
	found=true;
	continue;
      }
     }
     
     if(!found) {
      o->objects.push_back(*it); 
     }
   }
   return o;   
  }


  ///////////////////////////////////// BOOTSTRAP \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ObjectCollection* ObjectCollection::bootstrapWeight(string weightkey, double weightsum) const 
  {
   return bootstrapWeight(prototype->doubleVkey(weightkey),weightsum);
  }

  ObjectCollection* ObjectCollection::bootstrapWeight(int weightkey, double weightsum) const 
  {
       if(weightkey<0) {
  	cerr << "weighting with invalid double key" << endl;
        return 0;
       }

       ObjectCollection* o = new ObjectCollection();
       o->prototype = prototype; // same prototype

       double w = 0.;

       while(w<weightsum) {
         Object* obj = objects[rand()%size()];
         w += obj->doublePropertyValue(weightkey);
         o->objects.push_back(obj);
       }

       return o;
  }


  /////////////////////////////////// COLUMN OPERATIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      
  int ObjectCollection::createIntPropertyIfNecessary(string colname, int init, bool initevenifexists) // returns intVkey
  {
    int key=prototype->intVkey(colname);
    if(key>=0) {
      if(initevenifexists)
        setIntProperty(colname, init);
      return key;
    }
    prototype->intPropertyName.push_back(colname);
    prototype->intPropertyKey.push_back(-1); // no physical column key
    for(vector<Object*>::iterator it = objects.begin(); it != objects.end(); ++it) {
	(*it)->intProperty.push_back(init);
    }
    return prototype->intVkey(colname);
  }
  
  int ObjectCollection::createDoublePropertyIfNecessary(string colname, double init, bool initevenifexists) // returns doubleVkey
  {
    int key=prototype->doubleVkey(colname);
    if(key>=0) {
      if(initevenifexists)
        setDoubleProperty(colname, init);
      return key;
    }
    prototype->doublePropertyName.push_back(colname);
    prototype->doublePropertyKey.push_back(-1); // no physical column key
    for(vector<Object*>::iterator it = objects.begin(); it != objects.end(); ++it) {
	(*it)->doubleProperty.push_back(init);
    }
    return prototype->doubleVkey(colname);
  }  
  
  int ObjectCollection::createStringPropertyIfNecessary(string colname) // returns stringVkey
  {
    int key=prototype->stringVkey(colname);
    if(key>=0) {
      // no initialization
      return key;
    }
    prototype->stringPropertyName.push_back(colname);
    prototype->stringPropertyKey.push_back(-1); // no physical column key
    for(vector<Object*>::iterator it = objects.begin(); it != objects.end(); ++it) {
	(*it)->stringProperty.push_back("");
    }
    return prototype->stringVkey(colname);
  }
  
  bool ObjectCollection::calculateRp(string xcol, string ycol, double x0, double y0, string rpcol, string phipcol, bool clip_margin, double z0)
  {
    int rpkey=createDoublePropertyIfNecessary(rpcol);
    int phipkey=createDoublePropertyIfNecessary(phipcol);
    int xkey=prototype->doubleVkey(xcol);
    int ykey=prototype->doubleVkey(ycol);
    if(xkey<0 || ykey<0)
    {
      cerr << "x or y column not found, cannot calculate rp" << endl;
      return 0;
    }
   
    double rmax;
    if(clip_margin) { 
    	const double coord_min=0.;
    	const double coord_max=1.e6;
        rmax = sqrt(pow2_filter(min(coord_max-x0,x0-coord_min)) + pow2_filter(min(coord_max-y0,y0-coord_min)) + pow2_filter(min(coord_max-z0,z0-coord_min)));
    }

#pragma omp parallel for
    for(int i=0; i<objects.size(); i++) {
	objects[i]->doubleProperty[rpkey]=sqrt(pow2_filter(objects[i]->doubleProperty[xkey]-x0)+pow2_filter(objects[i]->doubleProperty[ykey]-y0));
	objects[i]->doubleProperty[phipkey]=atan2(objects[i]->doubleProperty[ykey]-y0,objects[i]->doubleProperty[xkey]-x0)+M_PI; // should be between 0 and 2*PI now
	if(clip_margin && objects[i]->doubleProperty[rpkey]>rmax) objects[i]->doubleProperty[rpkey]=-1.; // invalid projected radius
    }
    return 1;
  }
  
  bool ObjectCollection::calculateRpRADec(string racol, string deccol, double ra0, double dec0, string rpcol, string phipcol, double kpcperdeg)
  {
    int rpkey=createDoublePropertyIfNecessary(rpcol);
    int phipkey=createDoublePropertyIfNecessary(phipcol);
    int rakey=prototype->doubleVkey(racol);
    int deckey=prototype->doubleVkey(deccol);
    if(rakey<0 || deckey<0)
    {
      cerr << "RA or dec column not found, cannot calculate rp" << endl;
      return 0;
    }
    double c = cos(dec0*M_PI/180.);
#pragma omp parallel for
    for(int i=0; i<objects.size(); i++) {
	objects[i]->doubleProperty[rpkey]=kpcperdeg*sqrt(pow2_filter(c*(objects[i]->doubleProperty[rakey]-ra0))+pow2_filter(objects[i]->doubleProperty[deckey]-dec0));
	objects[i]->doubleProperty[phipkey]=atan2(dec0-objects[i]->doubleProperty[deckey],(ra0-objects[i]->doubleProperty[rakey])*c)+M_PI; // should be between 0 and 2*PI now
    }
    return 1;
  }
  
  bool ObjectCollection::calculateGtGx(string phicol, string g1col, string g2col, string gtcol, string gxcol)
  {
    int gtkey=createDoublePropertyIfNecessary(gtcol);
    int gxkey=createDoublePropertyIfNecessary(gxcol);
    int phikey = prototype->doubleVkey(phicol);
    int g1key = prototype->doubleVkey(g1col);
    int g2key = prototype->doubleVkey(g2col);
    if(phikey<0 || g1key<0 || g2key<0)
    {
      cerr << "Ellipticity or angle column not found, cannot calculate gt and gx" << endl;
      return 0;
    }
#pragma omp parallel for
    for(int i=0; i<objects.size(); i++) {
      double phi=objects[i]->doubleProperty[phikey];
      double c=cos(2.*phi); double s=sin(2.*phi);
      double g1=objects[i]->doubleProperty[g1key];
      double g2=objects[i]->doubleProperty[g2key];

      objects[i]->doubleProperty[gtkey]=-g1*c-g2*s;
      objects[i]->doubleProperty[gxkey]=-g1*s+g2*c;

//      if(i==0)
//      	cerr << "phi=" << phi << ", g1=" << g1 << ", g2=" << g2 << ", gt=" << objects[i]->doubleProperty[gtkey] << ", gx=" << objects[i]->doubleProperty[gxkey] << endl;

    }
    return 1;
  }
  
  double ObjectCollection::sum(string colname) const
  {
    double sum=0.;
    int ckey=prototype->doubleVkey(colname);
    if(ckey>=0) {
      for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	sum += (*it)->doubleProperty[ckey];
      }
      return sum;
    }
    else if((ckey=prototype->intVkey(colname)>=0)) {
      for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	sum += double((*it)->intProperty[ckey]);
      }
      return sum;
    }
    cerr << "# you are asking me to sum a non-existing column " << colname << " and I am rightfully refusing to" << endl;
    return -1.e99;
  }  
  double ObjectCollection::sumSq(string colname) const
  {
    double sum=0.;
    int ckey=prototype->doubleVkey(colname);
    if(ckey>=0) {
      for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	sum += pow((*it)->doubleProperty[ckey],2);
      }
      return sum;
    }
    else if((ckey=prototype->intVkey(colname)>=0)) {
      for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	sum += double(pow((*it)->intProperty[ckey],2));
      }
      return sum;
    }
    cerr << "# you are asking me to sum a non-existing column " << colname << " and I am rightfully refusing to" << endl;
    return -1.e99;
  }
  
  double ObjectCollection::sumProduct(string colname1, string colname2) const
  {
    double sum=0.;
    int ckey1=prototype->doubleVkey(colname1);
    int ckey2=prototype->doubleVkey(colname2);
    if(ckey1>=0 && ckey2>=0) {
      
#pragma omp parallel for reduction(+ : sum)
    for(int it=0; it<objects.size(); it++) {
    //  for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	double k = objects[it]->doubleProperty[ckey1]*objects[it]->doubleProperty[ckey2];
	sum += k;
      }
      return sum;
    }
    cerr << "double column not found in sumProduct" << endl;
    return 0.;
  }
  
  double ObjectCollection::sumErrSq(string col1, string col2) const
  {
    double ssq=0.;
    int ckey1=prototype->doubleVkey(col1);
    int ckey2=prototype->doubleVkey(col2);
    if(ckey1<0 || ckey2<0)
    {
     cerr << "sumErrSq: column " << col1 << " or " << col2 << " not found\nprototype: " << *prototype << endl; return 0;
    }
    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
      ssq += pow((*it)->doubleProperty[ckey1]-(*it)->doubleProperty[ckey2],2);
    }
    return ssq;
  }
  
  double ObjectCollection::sumErrSqWt(string col1, string col2, string col3) const
  {
    double ssq=0.;
    double swt=0.;
    int ckey1=prototype->doubleVkey(col1);
    int ckey2=prototype->doubleVkey(col2);
    int ckey3=prototype->doubleVkey(col3);
    if(ckey1<0 || ckey2<0 || ckey3<0)
    {
     cerr << "sumErrSq: column " << col1 << " or " << col2 << " or " << col3 << " not found\nprototype: " << *prototype << endl; return 0;
    }
    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
      ssq += pow((*it)->doubleProperty[ckey1]-(*it)->doubleProperty[ckey2],2)*(*it)->doubleProperty[ckey3];
      swt += (*it)->doubleProperty[ckey3];
    }
    return ssq/swt;
  }

  double ObjectCollection::sumErrSqWt(string col1, double col2, string col3) const
  {
    double ssq=0.;
    double swt=0.;
    int ckey1=prototype->doubleVkey(col1);
    int ckey3=prototype->doubleVkey(col3);
    if(ckey1<0 || ckey3<0)
    {
     cerr << "sumErrSq: column " << col1 << " or " << col3 << " not found\nprototype: " << *prototype << endl; return 0;
    }
    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
      ssq += pow((*it)->doubleProperty[ckey1]-col2,2)*(*it)->doubleProperty[ckey3];
      swt += (*it)->doubleProperty[ckey3];
    }
    return ssq/swt;
  }

  double ObjectCollection::sumErrSqIntObs(string col1, string col2, double interr, string obserr) const
  {
    double ssq=0.;
    int ckey1=prototype->doubleVkey(col1);
    int ckey2=prototype->doubleVkey(col2);
    int obserrkey=prototype->doubleVkey(obserr);
    if(ckey1<0 || ckey2<0 || obserrkey<0)
    {
     cerr << "sumErrSqIntObs: column " << col1 << " or " << col2 << " or " << obserr << " not found" << endl; return 0;
    }
    double interr2=interr*interr;
    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
      ssq += pow((*it)->doubleProperty[ckey1]-(*it)->doubleProperty[ckey2],2)/(interr2+pow((*it)->doubleProperty[obserrkey],2));
    }    
    return ssq;
  }
  
  double ObjectCollection::average(string colname) const
  {
    double sumt=sum(colname);
    return sumt/size();
  }
      
  double ObjectCollection::avgstdv(string colname) const
  {
    double sumt=sum(colname);
    double sum2t=sumSq(colname);
    double n = size();
    sumt /= n;
    sum2t /= n;
    return sqrt(sum2t-sumt*sumt)/sqrt(n-1.);
  }
    
  double ObjectCollection::averageWeighted(string colname, string weightname) const
  {
    double sum=0.;
    double wsum=0.;
    int ckey=prototype->doubleVkey(colname);
    int wkey=prototype->doubleVkey(weightname);
    if(ckey>=0 && wkey>=0) {
      for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	sum += (*it)->doubleProperty[ckey]*(*it)->doubleProperty[wkey];
	wsum += (*it)->doubleProperty[wkey];
      }
      return sum/wsum;
    }
    else if((ckey=prototype->intVkey(colname)>=0) && wkey>=0) {
      for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
	sum += double((*it)->intProperty[ckey])*(*it)->doubleProperty[wkey];
	wsum += (*it)->doubleProperty[wkey];
      }
      return sum/wsum;
    }
    return -1.e99;
  }
  
  double ObjectCollection::averageEllipsoid(string colname, vector<string> ecolnames, vector<double> centers, 
					    vector<double> radii, const vector<double> maxima, bool talk) const
  {
    //talk=true;
    
    vector<int> vkeys;
    vector<bool> outsidelimits;
    
    //cerr << "# average around " << centers[0] << " " << centers[1] << " " << centers[2] << endl;
    
    for(int i=0; i<ecolnames.size(); i++) {
     vkeys.push_back(prototype->doubleVkey(ecolnames[i]));
     if(vkeys[i]<0) {
      cerr << "invalid double property requested: " << ecolnames[i] << endl;
      return 0.;
     }
     if(maxima.size()==centers.size() && centers[i]>maxima[i]) outsidelimits.push_back(true);
     else outsidelimits.push_back(false);
    }
    
    int ckey=prototype->doubleVkey(colname);
    if(ckey<0) {
     cerr << "invalid double property requested: " << colname << endl;
    }
    
    double rfac=1.0;
    
    double sum;
    double wsum;
    int n;
    double a[ecolnames.size()];
    
    aeagain:
    
    sum=wsum=n=0.;
    for(int i=0; i<ecolnames.size(); i++) a[i]=0.;

    for(int it=0; it<objects.size(); it++) {
        double r2=0.;
	for(int i=0; i<ecolnames.size(); i++) {
	 double itmag=objects[it]->doubleProperty[vkeys[i]];
	 if(outsidelimits[i] && itmag>maxima[i]) {
	   // both tested center and reference object are above limits in this case, so nothing to add to the radius 
	 } else {
	   // tested center or reference object inside limits
	   r2 += pow2_filter((centers[i]-itmag)/radii[i]/rfac);
	 }
	 if(r2>=1) break;
	}
      
        if(r2<1) { // reference object is inside sphere

	  sum += objects[it]->doubleProperty[ckey]*(1.-r2);
	  for(int i=0; i<ecolnames.size(); i++) {
	    a[i] += objects[it]->doubleProperty[vkeys[i]]*(1.-r2);
	  }
	  wsum += (1.-r2);
	  //cerr << "adding " << objects[it]->doubleProperty[vkeys[0]] << " with weight " << 1.-r2 << endl;
	  //cerr << "average now " << a[0] << "/" << wsum << "=" << a[0]/wsum << endl;

	  n++;
        
	}
    }
    
    while(n<10 && rfac<3) {
     rfac *= 1.5;
     if(talk) cerr << "# too few objects, enlarging sphere by " << rfac << endl;
     goto aeagain;
    }
    
    if(n>=10)
    {
      // check for symmetry of sample
      for(int i=0; i<ecolnames.size(); i++) {
	if(!outsidelimits[i] && fabs(a[i]/wsum-centers[i])>2.*radii[i]) { // only if inside limits, check for shift
	  cerr << "# " << a[i]/wsum << " is far from " << centers[i] << " in terms of " << ecolnames[i] << endl;
	  return 0./0.;
	}
      }
      return sum/wsum;
    }

    return 0./0.;
  }
        
  Object* ObjectCollection::randomFromEllipsoid(vector<string> ecolnames, vector<double> centers, vector<double> radii) const
  {
    vector<int> vkeys;
    for(int i=0; i<ecolnames.size(); i++) {
     vkeys.push_back(prototype->doubleVkey(ecolnames[i]));
     if(vkeys[i]<0) {
      cerr << "invalid double property requested: " << ecolnames[i] << endl;
      return 0;
     }
    }
    
    double rfac=1.0;
    
    veagain:
    double sum=0.;
    double wsum=0.;

    vector<Object *> vval;
    
#pragma omp parallel for
    for(int it=0; it<objects.size(); it++) {
//    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
        double r2=0.;
	for(int i=0; i<ecolnames.size(); i++) {
	 r2 += pow2_filter((centers[i]-objects[it]->doubleProperty[vkeys[i]])/radii[i]/rfac);
	 if(r2>=1) break;
	}
      
        if(r2<1) {
#pragma omp critical
	  vval.push_back(objects[it]);
	}
    }
    
    while(vval.size()<1 && rfac<4) {
     rfac *= 1.2;
     goto veagain;
    }
    if(vval.size()==0) return 0;
    
    return vval[rand()%vval.size()];
  }
  
  int ObjectCollection::countInEllipsoid(vector<string> ecolnames, vector<double> centers, vector<double> radii) const
  {
    vector<int> vkeys;
    for(int i=0; i<ecolnames.size(); i++) {
     vkeys.push_back(prototype->doubleVkey(ecolnames[i]));
     if(vkeys[i]<0) {
      cerr << "invalid double property requested: " << ecolnames[i] << endl;
      return 0;
     }
    }
  
    int vc=0;
    
#pragma omp parallel for
    for(int it=0; it<objects.size(); it++) {
//    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
        double r2=0.;
	for(int i=0; i<ecolnames.size(); i++) {
	 r2 += pow2_filter((centers[i]-objects[it]->doubleProperty[vkeys[i]])/radii[i]);
	 if(r2>=1) break;
	}
      
        if(r2<1) {
#pragma omp critical
	  vc++;
	}
    }
    
    return vc;
  }
  
   
  double ObjectCollection::percentile(string colname, double p)
  {
    int dckey=prototype->doubleVkey(colname);
    int ickey=prototype->intVkey(colname);
    if(dckey<0 && ickey<0)
    {
     cerr << "ERROR: cannot determine median of unknown double or int column " << colname << endl;
     return 0./0.;
    }
    int n = objects.size();
    //cerr << "# determining median of " << colname << " from " << n << " objects" << endl; 
    sort(colname);

    return objects[int(p*double(n))+1]->numericValue(colname);
    
  }
  
  double ObjectCollection::median(string colname)
  {
    int dckey=prototype->doubleVkey(colname);
    int ickey=prototype->intVkey(colname);
    if(dckey<0 && ickey<0)
    {
     cerr << "ERROR: cannot determine median of unknown double or int column " << colname << endl;
     return 0./0.;
    }
    int n = objects.size();
    //cerr << "# determining median of " << colname << " from " << n << " objects" << endl; 
    sort(colname);
    
    
    if(n%2==1) // odd number
    {
     return objects[(n-1)/2]->numericValue(colname);
    }
    return (objects[n/2]->numericValue(colname) + objects[n/2-1]->numericValue(colname))/2.;
  }

  void ObjectCollection::transformColumn(string colname, double (*f)(double)) {
    int ckey=prototype->doubleVkey(colname);
    if(ckey<0) {
     ckey=prototype->intVkey(colname);
     if(ckey<0) {
      cerr << "ERROR: cannot transform unknown double or int column " << colname << endl;
      return;
     } else { // found an int, not a double
#pragma omp parallel for
	for(int iit=0; iit<objects.size(); iit++) {
	Object **it= &(objects[iit]);
	(*it)->intProperty[ckey]=(*f)((*it)->intProperty[ckey]);
	}
	return;
     }
    }
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
     (*it)->doubleProperty[ckey]=(*f)((*it)->doubleProperty[ckey]);
    }
  }
   
  void ObjectCollection::transformColumn(string colname, double (*f)(double, double), string arg1) {
    int ckey=prototype->doubleVkey(colname);
    int ckey_int=prototype->intVkey(colname);
    int a1=prototype->doubleVkey(arg1);
    int a1_int=prototype->intVkey(arg1);
    if((ckey<0 && ckey_int<0) || (a1<0 && a1_int<0)) {
     cerr << "ERROR: cannot transform unknown column" << endl;
     if(ckey<0 && ckey_int<0) cerr << "missing: " << colname << endl;
     if(a1<0 && a1_int<0) cerr << "missing: " << arg1 << endl;
     cerr << "prototype: " << *prototype << endl;
     return;
    }
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
    if(ckey>=0) { // transformed key is double
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
     (*it)->doubleProperty[ckey]=(*f)((*it)->doubleProperty[ckey],(a1>=0)?((*it)->doubleProperty[a1]):(double((*it)->intProperty[a1_int])));
    }
    } else { // transformed property is int
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
     (*it)->intProperty[ckey_int]=int((*f)(double((*it)->intProperty[ckey_int]),(a1>=0)?((*it)->doubleProperty[a1]):(double((*it)->intProperty[a1_int]))));
    }
    }
  }
     
  void ObjectCollection::transformColumn(string colname, double (*f)(double, double, double), string arg1, string arg2) {
    int ckey=prototype->doubleVkey(colname);
    int ckey_int=prototype->intVkey(colname);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    if((ckey<0 && ckey_int<0) || a1<0 || a2<0) {
     cerr << "ERROR: cannot transform unknown double column" << endl;
     if(ckey<0 && ckey_int<0) cerr << "missing: " << colname << endl;
     if(a1<0) cerr << "missing: " << arg1 << endl;
     if(a2<0) cerr << "missing: " << arg2 << endl;
     cerr << "prototype: " << *prototype << endl;
     return;
    }
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
    if(ckey>=0) { // transformed key is double
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
     (*it)->doubleProperty[ckey]=(*f)((*it)->doubleProperty[ckey],(*it)->doubleProperty[a1],(*it)->doubleProperty[a2]);
    }
    } else { // transformed property is int
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
     (*it)->intProperty[ckey_int]=int((*f)(double((*it)->intProperty[ckey_int]),(*it)->doubleProperty[a1], (*it)->doubleProperty[a2]));
    }
    }
  }     
  
  void ObjectCollection::transformColumn(string colname, double (*f)(double, double, double, double), string arg1, string arg2, string arg3) {
    int ckey=prototype->doubleVkey(colname);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    int a3=prototype->doubleVkey(arg3);
    if(ckey<0 || a1<0 || a2<0 || a3<0) {
     cerr << "ERROR: cannot transform unknown double column" << endl;
     if(ckey<0) cerr << "missing: " << colname << endl;
     if(a1<0) cerr << "missing: " << arg1 << endl;
     if(a2<0) cerr << "missing: " << arg2 << endl;
     if(a3<0) cerr << "missing: " << arg3 << endl;
     cerr << "prototype: " << *prototype << endl;
     return;
    }
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
     (*it)->doubleProperty[ckey]=(*f)((*it)->doubleProperty[ckey],(*it)->doubleProperty[a1],(*it)->doubleProperty[a2],(*it)->doubleProperty[a3]);
    }
  }
  
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double), string arg) {
    int ckey=createDoublePropertyIfNecessary(newcolname, -1, false);
    int a1=prototype->doubleVkey(arg);
    int a1_int=prototype->intVkey(arg);
    if(a1<0 && a1_int<0) {
      cerr << "ERROR: cannot use unknown double or int column for transformation" << endl;
      if(a1<0 && a1_int<0) cerr << "missing: " << arg << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
    
    if(a1>=0) { // found a double column
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1]);
    }
    } else { // found an int column
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->intProperty[a1_int]);
    }      
    }
    return ckey;
  }
  
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double), string arg1, string arg2) {
    int ckey=createDoublePropertyIfNecessary(newcolname, -1, false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    if(a1<0 || a2<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2]);
    }
    return ckey;
  }
    
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double), string arg1, string arg2, string arg3) {
    int ckey=createDoublePropertyIfNecessary(newcolname, -1, false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    int a3=prototype->doubleVkey(arg3);
    if(a1<0 || a2<0 || a3<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3]);
    }
    return ckey;
  }
  
    
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double), string arg1, string arg2, string arg3, string arg4) {
    int ckey=createDoublePropertyIfNecessary(newcolname, -1, false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    if(a1<0 || a2<0 || a3<0 || a4<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4]);
    }
    return ckey;
  }
  
      
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5) {
    int ckey=createDoublePropertyIfNecessary(newcolname, -1, false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5]);
    }
    return ckey;
  }
  
        
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6]);
    }
    return ckey;
  }
          
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2);
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7]);
    }
    return ckey;
  }

            
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8]);
    }
    return ckey;
  }
  
              
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9]);
    }
    return ckey;
  }
  
              
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10]);
    }
    return ckey;
  }

                
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    int a11=prototype->doubleVkey(arg11);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0 || a11<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      if(a11<0) cerr << "missing: " << arg11 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10],(*it)->doubleProperty[a11]);
    }
    return ckey;
  }
                
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    int a11=prototype->doubleVkey(arg11);
    int a12=prototype->doubleVkey(arg12);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0 || a11<0 || a12<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      if(a11<0) cerr << "missing: " << arg11 << endl;
      if(a12<0) cerr << "missing: " << arg12 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10],(*it)->doubleProperty[a11],(*it)->doubleProperty[a12]);
    }
    return ckey;
  }
  
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double,double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    int a11=prototype->doubleVkey(arg11);
    int a12=prototype->doubleVkey(arg12);
    int a13=prototype->doubleVkey(arg13);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0 || a11<0 || a12<0 || a13<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      if(a11<0) cerr << "missing: " << arg11 << endl;
      if(a12<0) cerr << "missing: " << arg12 << endl;
      if(a13<0) cerr << "missing: " << arg13 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10],(*it)->doubleProperty[a11],(*it)->doubleProperty[a12],(*it)->doubleProperty[a13]);
    }
    return ckey;
  }
 
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double,double,double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13, string arg14) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    int a11=prototype->doubleVkey(arg11);
    int a12=prototype->doubleVkey(arg12);
    int a13=prototype->doubleVkey(arg13);
    int a14=prototype->doubleVkey(arg14);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0 || a11<0 || a12<0 || a13<0 || a14<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      if(a11<0) cerr << "missing: " << arg11 << endl;
      if(a12<0) cerr << "missing: " << arg12 << endl;
      if(a13<0) cerr << "missing: " << arg13 << endl;
      if(a14<0) cerr << "missing: " << arg14 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10],(*it)->doubleProperty[a11],(*it)->doubleProperty[a12],(*it)->doubleProperty[a13],(*it)->doubleProperty[a14]);
    }
    return ckey;
  }
  
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double,double,double,double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13, string arg14, string arg15) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    int a11=prototype->doubleVkey(arg11);
    int a12=prototype->doubleVkey(arg12);
    int a13=prototype->doubleVkey(arg13);
    int a14=prototype->doubleVkey(arg14);
    int a15=prototype->doubleVkey(arg15);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0 || a11<0 || a12<0 || a13<0 || a14<0 || a15<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      if(a11<0) cerr << "missing: " << arg11 << endl;
      if(a12<0) cerr << "missing: " << arg12 << endl;
      if(a13<0) cerr << "missing: " << arg13 << endl;
      if(a14<0) cerr << "missing: " << arg14 << endl;
      if(a15<0) cerr << "missing: " << arg15 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10],(*it)->doubleProperty[a11],(*it)->doubleProperty[a12],(*it)->doubleProperty[a13],(*it)->doubleProperty[a14],(*it)->doubleProperty[a15]);
    }
    return ckey;
  }
   
  int ObjectCollection::transformColumnNew(string newcolname, double (*f)(double,double,double,double,double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13, string arg14, string arg15, string arg16) {
    int ckey=createDoublePropertyIfNecessary(newcolname,-1,false);
    int a1=prototype->doubleVkey(arg1);
    int a2=prototype->doubleVkey(arg2); 
    int a3=prototype->doubleVkey(arg3);
    int a4=prototype->doubleVkey(arg4);
    int a5=prototype->doubleVkey(arg5);
    int a6=prototype->doubleVkey(arg6);
    int a7=prototype->doubleVkey(arg7);
    int a8=prototype->doubleVkey(arg8);
    int a9=prototype->doubleVkey(arg9);
    int a10=prototype->doubleVkey(arg10);
    int a11=prototype->doubleVkey(arg11);
    int a12=prototype->doubleVkey(arg12);
    int a13=prototype->doubleVkey(arg13);
    int a14=prototype->doubleVkey(arg14);
    int a15=prototype->doubleVkey(arg15);
    int a16=prototype->doubleVkey(arg16);
    if(a1<0 || a2<0 || a3<0 || a4<0 || a5<0 || a6<0 || a7<0 || a8<0 || a9<0 || a10<0 || a11<0 || a12<0 || a13<0 || a14<0 || a15<0 || a16<0) {
      cerr << "ERROR: cannot use unknown double column for transformation" << endl;
      if(a1<0) cerr << "missing: " << arg1 << endl;
      if(a2<0) cerr << "missing: " << arg2 << endl;
      if(a3<0) cerr << "missing: " << arg3 << endl;
      if(a4<0) cerr << "missing: " << arg4 << endl;
      if(a5<0) cerr << "missing: " << arg5 << endl;
      if(a6<0) cerr << "missing: " << arg6 << endl;
      if(a7<0) cerr << "missing: " << arg7 << endl;
      if(a8<0) cerr << "missing: " << arg8 << endl;
      if(a9<0) cerr << "missing: " << arg9 << endl;
      if(a10<0) cerr << "missing: " << arg10 << endl;
      if(a11<0) cerr << "missing: " << arg11 << endl;
      if(a12<0) cerr << "missing: " << arg12 << endl;
      if(a13<0) cerr << "missing: " << arg13 << endl;
      if(a14<0) cerr << "missing: " << arg14 << endl;
      if(a15<0) cerr << "missing: " << arg15 << endl;
      if(a16<0) cerr << "missing: " << arg16 << endl;
      cerr << "prototype: " << *prototype << endl;
      return -1;
    }
     
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
#pragma omp parallel for
    for(int iit=0; iit<objects.size(); iit++) {
     Object **it= &(objects[iit]);
      (*it)->doubleProperty[ckey] = (*f)((*it)->doubleProperty[a1], (*it)->doubleProperty[a2], (*it)->doubleProperty[a3], (*it)->doubleProperty[a4], (*it)->doubleProperty[a5], (*it)->doubleProperty[a6], (*it)->doubleProperty[a7], (*it)->doubleProperty[a8],(*it)->doubleProperty[a9],(*it)->doubleProperty[a10],(*it)->doubleProperty[a11],(*it)->doubleProperty[a12],(*it)->doubleProperty[a13],(*it)->doubleProperty[a14],(*it)->doubleProperty[a15],(*it)->doubleProperty[a16]);
    }
    return ckey;
  }
  
  void ObjectCollection::setIntProperty(string colname, int value, int first, int last)
  {
     int ckey=prototype->intVkey(colname);
     if(ckey<0) {
      cerr << "ERROR: cannot set unknown int column " << colname << endl;
      return; 
     }
     for(int i=first; i<((last>=0 && last<objects.size())?(last+1):objects.size()); i++)
     {
      objects[i]->intProperty[ckey]=value; 
     }
  }
  
  void ObjectCollection::setDoubleProperty(string colname, double value, int first, int last)
  {
     int ckey=prototype->doubleVkey(colname);
     if(ckey<0) {
      cerr << "ERROR: cannot set unknown double column " << colname << endl;
      return; 
     }
     //cerr << "setting double property " << colname << " to " << value << endl;
     for(int i=first; i<((last>=0 && last<objects.size())?(last+1):objects.size()); i++)
     {
      objects[i]->doubleProperty[ckey]=value; 
     }
  }
    
  ///////////////////////////////////////// BINNING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
void ObjectCollection::sort(string sortcolumn)
{
    prototype->setSortKey(sortcolumn);
    std::sort(objects.begin(), objects.end(), ObjectPrototype::sortObjects);
}
    
ObjectCollection* ObjectCollection::bin(vector<string> averaged_columns, string binning_column, vector<double> binning_limits)
{
  // sort by binning key
  sort(binning_column);
  
  // prepare binned table
  ObjectCollection *bins = new ObjectCollection();
  string *binnames = new string[2*averaged_columns.size()+2];
  int *binkeys = new int[2*averaged_columns.size()+2];
  Type *bintypes = new Type[2*averaged_columns.size()+2];
  for(int i=0; i<2*averaged_columns.size()+2; i++)
  {
      binkeys[i]=-1;
      bintypes[i]=doubleType;
  }
  
  binnames[0] = binning_column+"_min";
  binnames[1] = binning_column+"_max";
  int okey[averaged_columns.size()+1]; // vkeys of the columns to be averaged
  Type otype[averaged_columns.size()+1]; // types of the columns to be averaged
  okey[0]=prototype->intVkey(binning_column);
  otype[0]=integerType;
  if(okey[0]<0)
  {
    okey[0]=prototype->doubleVkey(binning_column);
    otype[0]=doubleType;
  }
  if(okey[0]<0)
  {
    cerr << "ERROR: sorting type is neither integer nor double!" << endl;
    return 0;
  }
  
  for (int i=2; i<averaged_columns.size()+2; i++)
  {
   binnames[i]=averaged_columns[i-2];
   binnames[i+averaged_columns.size()]="sig_"+averaged_columns[i-2];
   if((okey[i-1]=prototype->intVkey(averaged_columns[i-2]))>-1)
    otype[i-1]=integerType;
   else if((okey[i-1]=prototype->doubleVkey(averaged_columns[i-2]))>-1)
    otype[i-1]=doubleType;
   else
   {
    cerr << "ERROR: averaged type is neither integer nor double!" << endl;
    return 0;
   }
  }
  bins->setPrototype(binkeys, binnames, bintypes, 2*averaged_columns.size()+2);
  bins->makeEmpty(binning_limits.size()-1);
  
  // average the properties in bins
  int omin=-1;
  {
  int i=0;
  while(omin<0 && i<objects.size()) {
	if(binning_limits[0]<=objects[i]->doubleProperty[okey[0]])
	{
	  omin=i;
	}
	i++;
  }
  }
  int omax=omin;
  
  for (int j=0; j<binning_limits.size()-1; j++)
  {
    {
    int i=omin;
    while(omax==omin && i<objects.size()) {
      	if(binning_limits[j+1]<=objects[i]->doubleProperty[okey[0]])
	{
	  omax=i-1;
	}
	i++;
    }
    }

    (*bins)[j]->doubleProperty[0] = binning_limits[j];
    (*bins)[j]->doubleProperty[1] = binning_limits[j+1];

    //cerr << "averaging between id " << omin << " and " << omax << endl;
    
    for (int i=2; i<averaged_columns.size()+2; i++)
    {
	double sum=0.;
	double sum2=0.;
	if(otype[i-1]==integerType) {
	 for (int k=omin; k<=omax; k++) {
	   int o = objects[k]->intProperty[okey[i-1]];
	   sum += double(o);
	   sum2+= double(o*o);
	 }
	}
	else if (otype[i-1]==doubleType) {
	 for (int k=omin; k<=omax; k++)
	 {
	  double o = objects[k]->doubleProperty[okey[i-1]];
	  sum  += o;
	  sum2 += o*o;
	 }
	}
	(*bins)[j]->doubleProperty[i] = sum / double(omax-omin+1);
	(*bins)[j]->doubleProperty[averaged_columns.size()+i] 
		    = sqrt(1./double(omax-omin)*((sum2 / double(omax-omin+1))-pow2_filter((*bins)[j]->doubleProperty[i])));
    }
    
    omin=omax+1;
    omax=omin;
  }
  
  return bins;
}


    
ObjectCollection* ObjectCollection::bin(vector<string> averaged_columns, string binning_column, int nbins)
    // average averaged_columns for nbins bins in binning_column, of equal population 
    // warning: may change order of this ObjectCollection; the returned table is allocated and must be disallocated by the user
{
  // sort by binning key
  sort(binning_column);
  
  // prepare binned table
  ObjectCollection *bins = new ObjectCollection();
  string *binnames = new string[2*averaged_columns.size()+2];
  int *binkeys = new int[2*averaged_columns.size()+2];
  Type *bintypes = new Type[2*averaged_columns.size()+2];
  for(int i=0; i<2*averaged_columns.size()+2; i++)
  {
      binkeys[i]=-1;
      bintypes[i]=doubleType;
  }
  
  binnames[0] = binning_column+"_min";
  binnames[1] = binning_column+"_max";
  int okey[averaged_columns.size()+1]; // vkeys of the columns to be averaged
  Type otype[averaged_columns.size()+1]; // types of the columns to be averaged
  okey[0]=prototype->intVkey(binning_column);
  otype[0]=integerType;
  if(okey[0]<0)
  {
    okey[0]=prototype->doubleVkey(binning_column);
    otype[0]=doubleType;
  }
  if(okey[0]<0)
  {
    cerr << "ERROR: sorting type is neither integer nor double!" << endl;
    return 0;
  }
  
  for (int i=2; i<averaged_columns.size()+2; i++)
  {
   binnames[i]=averaged_columns[i-2];
   binnames[i+averaged_columns.size()]="sig_"+averaged_columns[i-2];
   if((okey[i-1]=prototype->intVkey(averaged_columns[i-2]))>-1)
    otype[i-1]=integerType;
   else if((okey[i-1]=prototype->doubleVkey(averaged_columns[i-2]))>-1)
    otype[i-1]=doubleType;
   else
   {
    cerr << "ERROR: averaged type is neither integer nor double!" << endl;
    return 0;
   }
  }
  bins->setPrototype(binkeys, binnames, bintypes, 2*averaged_columns.size()+2);
  bins->makeEmpty(nbins);
  
  // average the properties in bins
  int omax=-1;
  for (int j=0; j<nbins; j++)
  {
    int omin=omax+1;				// first object in bin
    omax=objects.size()/nbins*(j+1)-1;		//  last object in bin
    if(j==nbins-1) omax=objects.size()-1;

    (*bins)[j]->doubleProperty[0] = objects[omin]->doubleProperty[okey[0]];
    (*bins)[j]->doubleProperty[1] = objects[omax]->doubleProperty[okey[0]];

    for (int i=2; i<averaged_columns.size()+2; i++)
    {
	double sum=0.;
	double sum2=0.;
	if(otype[i-1]==integerType) {
	 for (int k=omin; k<=omax; k++) {
	   int o = objects[k]->intProperty[okey[i-1]];
	   sum += double(o);
	   sum2+= double(o*o);
	 }
	}
	else if (otype[i-1]==doubleType) {
	 for (int k=omin; k<=omax; k++)
	 {
	  double o = objects[k]->doubleProperty[okey[i-1]];
	  sum  += o;
	  sum2 += o*o;
	 }
	}
	(*bins)[j]->doubleProperty[i] = sum / double(omax-omin+1);
	(*bins)[j]->doubleProperty[averaged_columns.size()+i] 
		    = sqrt(1./double(omax-omin)*((sum2 / double(omax-omin+1))-pow2_filter((*bins)[j]->doubleProperty[i])));
    }
  }
  
  return bins;
}

ObjectCollection* ObjectCollection::bin(vector<string> averaged_columns, vector<string> binning_columns, vector<int> nbins)
{
  int dim = binning_columns.size();
  if(nbins.size()!=dim) {
   cerr << "ERROR: did not pass a consistent binning dimensionality!" << endl;
   return 0;
  }
  
  ObjectCollection *bintable = bin(averaged_columns, binning_columns[0], nbins[0]);
  
  for(int n=1; n<dim; n++)
  {
    ObjectCollection *oldbintable = bintable;
    bintable = 0;
    // take each bin, filter by its limits, split again in new bins
    for(int i=0; i<oldbintable->size(); i++)
    {
     vector<Filter> selectbin;
     for(int j=0; j<n; j++)
     {
      selectbin.push_back(Filter(binning_columns[j],((*oldbintable)[i])->doublePropertyValue(binning_columns[j]+"_min"),((*oldbintable)[i])->doublePropertyValue(binning_columns[j]+"_max")));
     }
     ObjectCollection *subsample = filter(selectbin);
     ObjectCollection *subbintable = subsample->bin(averaged_columns, binning_columns[n], nbins[n]);
     delete subsample;
     for(int j=0; j<n; j++)
     {
      subbintable->createDoublePropertyIfNecessary(binning_columns[j]+"_min",((*oldbintable)[i])->doublePropertyValue(binning_columns[j]+"_min"));
      subbintable->createDoublePropertyIfNecessary(binning_columns[j]+"_max",((*oldbintable)[i])->doublePropertyValue(binning_columns[j]+"_max"));
     }
     if(!bintable) {
      bintable = subbintable; // don't delete subbintable
     } else {
      bintable->appendCollectionDeep(*subbintable); 
      delete subbintable;
     }
    }
    // bintable now contains all n+1-dimensional bins
    // oldbintable is not required any more
    delete oldbintable;
  }
  
  return bintable; 
}
    
void ObjectCollection::applyBinAverage(ObjectCollection *bintable, string column, vector<string> binning_columns, string binaveragecolumn_name)
{
  int tbkey=bintable->prototype->doubleVkey(column);
  int obkey=createDoublePropertyIfNecessary(binaveragecolumn_name);

  vector<int> binning_column_minvkeys;
  vector<int> binning_column_maxvkeys;
  vector<int> object_column_vkeys;
  for(vector<string>::const_iterator s = binning_columns.begin(); s != binning_columns.end(); ++s) {
    binning_column_minvkeys.push_back(bintable->prototype->doubleVkey(*s+"_min"));
    binning_column_maxvkeys.push_back(bintable->prototype->doubleVkey(*s+"_max"));
    object_column_vkeys.push_back(prototype->doubleVkey(*s));
  }
  
  for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
    int b = 0; // which bin are we in?
    int inbin = -1;
    while (inbin<0 && b<bintable->size()) {
      //cerr << "testing bin " << b+1 << " / " << bintable->size() << endl;
      for(int s = 0; s < binning_columns.size(); s++) {
	//cerr << "testing criterion " << s+1 << " / " << binning_columns.size() << endl;
	//cerr << (*bintable)[b]->doublePropertyValue(binning_column_minvkeys[s]) << " < " << (*it)->doublePropertyValue(object_column_vkeys[s]) << " < " << (*bintable)[b]->doublePropertyValue(binning_column_maxvkeys[s]) << endl;
	if(  (*bintable)[b]->doublePropertyValue(binning_column_minvkeys[s]) > (*it)->doublePropertyValue(object_column_vkeys[s]) ||
	     (*bintable)[b]->doublePropertyValue(binning_column_maxvkeys[s]) < (*it)->doublePropertyValue(object_column_vkeys[s])   ) {     
	     //cerr << "no!" << endl;
	     goto outofbin;
	}
	//cerr << "yes!" << endl;
      }
      //cerr << "accepting bin" << endl;
      inbin=b;
      outofbin:
      b++;
    } 
    
    if (inbin<0)
    {
     cerr << "ERROR: object does not fit in any bin!" << endl; 
     return;
    }
    //cerr << "setting binned property" << endl;
    (*it)->setDoubleProperty(obkey,(*bintable)[inbin]->doublePropertyValue(tbkey));
  }
}
  
  /////////////////////////////////// HELPFUL SMALL STUFF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
  int ObjectCollection::size() const 
  {
   return objects.size(); 
  }
  
  Object* ObjectCollection::operator[] (const int i) const
  {
   if(i>=0 && i<objects.size())
    return objects[i]; 
   return 0;
  }
  
  void ObjectCollection::print(vector<string> column, ostream& sink, bool header) const
  {
    if(header) {
      sink << "#";
    }
      vector<int> vkeys;
      vector<Type> types;
    for(vector<string>::iterator it = column.begin(); it != column.end(); it++) {
	if(header) sink << *it << " ";
	int i=prototype->doubleVkey(*it);
	if(i>=0) { types.push_back(doubleType); vkeys.push_back(i); }
	else { i=prototype->intVkey(*it); if(i<0) {cerr << "cannot print unavailable column " << *it << endl; return; } types.push_back(integerType); vkeys.push_back(i); }
    }
    if(header) sink << endl;
    sink.precision(10);
    sink.setf(ios::fixed,ios::floatfield);
    for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
        for(int s = 0; s < column.size(); s++) {
	  sink << ((types[s]==doubleType)?((*it)->doubleProperty[vkeys[s]]):((*it)->intProperty[vkeys[s]])) << " ";
	}
	sink << endl;
    }
  }
  
      
//////////////////////////////////////// SMART STUFF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
double ObjectCollection::epanechnikovDensity(string cx, string cy, double x, double y, double r) const { // 2d epanechnikov kernel density with radius r at x,y
    int ix = prototype->doubleVkey(cx);
    int iy = prototype->doubleVkey(cy);
    if(ix<0 || iy<0) {
     cerr << "x or y column not found" << endl; return 0.; 
    }
    
    double sum=0.;
    int n=0;

#pragma omp parallel for 
    //for(vector<Object*>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
    for(int i=0; i<objects.size(); i++) {
        Object *it = objects[i];
        double ro = sqrt(pow(x-it->doubleProperty[ix],2)+pow(y-it->doubleProperty[iy],2));
	if(ro<r) {
           double d = (1.-pow(ro/r,2));
#pragma omp critical
{
	   sum += d;
           n++;
}
        }
    }
    //cerr << "# " << n/3.14159 << " " << 2.*sum/3.14159 << endl;
    return 2.*sum/3.14159; // 2d normalization to objects per r^2: int( (1-r^2) 2pi r dr ) = pi/2 for unity density
}

ostream& operator<<(ostream& os, Object& dt)
{
    for(vector<int>::iterator it = dt.intProperty.begin(); it != dt.intProperty.end(); ++it) {
	os << *it << ' '; 
    }
    
    for(vector<double>::iterator it = dt.doubleProperty.begin(); it != dt.doubleProperty.end(); ++it) {
	os << *it << ' '; 
    }
    
    for(vector<string>::iterator it = dt.stringProperty.begin(); it != dt.stringProperty.end(); ++it) {
	os << *it << ' '; 
    }
    return os;
}

ostream& operator<<(ostream &os, ObjectPrototype& prototype)
{
    os << "# ";
  
    for(vector<string>::iterator it = prototype.intPropertyName.begin(); it != prototype.intPropertyName.end(); ++it) {
	os << *it << " "; 
    }
    
    for(vector<string>::iterator it = prototype.doublePropertyName.begin(); it != prototype.doublePropertyName.end(); ++it) {
	os << *it << " "; 
    }
    
    for(vector<string>::iterator it = prototype.stringPropertyName.begin(); it != prototype.stringPropertyName.end(); ++it) {
	os << *it << " "; 
    }
    
    os << '\n';
  
}

ostream& operator<<(ostream& os, ObjectCollection& dt)
{
    os << *(dt.prototype);
    
    for(vector<Object*>::iterator it = dt.objects.begin(); it != dt.objects.end(); ++it) {
	os << *(*it) << '\n'; 
    }
    
    return os;
}

