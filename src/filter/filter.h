/**** filter.h ****/
/**** Daniel Gruen, 2013-2014 ****/

/** These classes provide a multi-threaded library for catalog operations.
 *  The main class is an ObjectCollection, which contains a vector of Objects
 *  Objects are a vector of properties, ordered and labelled according to 
 *  an ObjectPrototype.
 *  An ObjectCollection contain methods for bulk operations on objects, such as
 *  filtering, matching, averaging, density estimation, the transformation of 
 *  known properties into new columns, the concatenation with other collections.
 *  They can be initiated from ASCII or FITS tables.
 *  Objects contain little more than access functions to their integer, double
 *  and string properties.
 * **/

#ifndef FILTER_H
#define FILTER_H

#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>

#ifndef NO_OMP
#include <omp.h>
#endif

using namespace std;


// class Filter:
//   a filter condition for some property to be between two values (inclusive >= and <=)
class Filter {
public:
 Filter(string _name, double _minVal, double _maxVal): name(_name), minVal(_minVal), maxVal(_maxVal)
 { }
 Filter(string _name, vector<double> _removeSet): name(_name), removeSet(_removeSet)
 { }
 string name;
 double minVal;
 double maxVal;
 vector<double> removeSet;
};

enum Type
{
  integerType,
  doubleType,
  stringType
};

class ObjectPrototype;

// class Object:
//   holds information about one object, including a pointer to its prototype and vectors of all properties
//   the best philosophy is to do as little as possible on an Object level - ObjectCollection has many useful 
//   bulk operations implemented in a parallelized fashion already

class Object {
  friend class ObjectCollection;
  friend class ObjectPrototype;
  public:
    
  Object(ObjectPrototype *p);			// makes an empty object with all properties set to 0
  Object(ObjectPrototype *p, ifstream &file);	// also reads object from text file
  
  int intPropertyValue(int vkey) const;		// return or set integer property based on its column key or name
  int intPropertyValue(string name) const;      
  bool setIntProperty(int vkey, int value);
  bool setIntProperty(string name, int value);
  
  double doublePropertyValue(int vkey) const;   // return or set double property based on its column key or name
  double doublePropertyValue(string name) const;
  bool setDoubleProperty(int vkey, double value);
  bool setDoubleProperty(string name, double value);
  
  string stringPropertyValue(int vkey) const;   // return or set string property based on its column key or name
  string stringPropertyValue(string name) const;
  bool setStringProperty(int vkey, string value);
  bool setStringProperty(string name, string value);  
  
  double numericValue(string name) const;       // a convenience wrapper for double or integer properties

  private:
      
  void readLine(ifstream &file);		// read object from line from file
  
  void makeEmpty();				// set all properties to 0 / empty string
  
  void readObject(string line);			// read object from line
  
  // Property vectors
  
  vector<int>    intProperty;
  vector<double> doubleProperty;
  vector<string> stringProperty;
  
  ObjectPrototype *prototype;
  
  static bool isComment(const string& instr) { // true if line is comment
    std::istringstream is(instr);
    string word1;
    is >> word1;
    return (!is || word1.empty() || word1[0]=='#');
  }

  static ifstream& getlineNoComment(ifstream& is, string& s) {
    do {
      if (!getline(is,s)) return is;
    } while (isComment(s));
    return is;
  }
  
  friend ostream& operator<<(ostream& os, Object& dt);
};

// class ObjectPrototype:
//   holds header information about what an object looks like,
//   in particular vectors of names and columns of all integer, double and string properties separately
//   you should not have to mess with this unless working on new methods in the Object and ObjectCollection classes

class ObjectPrototype {
    
  public:
    
    ObjectPrototype(): intSortVkey(-1), doubleSortVkey(-1) { }; // empty prototype
    ObjectPrototype(ObjectPrototype *pp): intSortVkey(-1), doubleSortVkey(-1), intPropertyName(pp->intPropertyName), doublePropertyName(pp->doublePropertyName), stringPropertyName(pp->stringPropertyName), intPropertyKey(pp->intPropertyKey), doublePropertyKey(pp->doublePropertyKey), stringPropertyKey(pp->stringPropertyKey) { };
    ObjectPrototype(const int keys[], const string names[], const Type types[], int n, string prefix="");
    
    vector<string> intPropertyName;
    vector<string> doublePropertyName;
    vector<string> stringPropertyName;
    vector<int> intPropertyKey;		// denotes column number for ascii input, no meaning otherwise
    vector<int> doublePropertyKey;
    vector<int> stringPropertyKey;
    
    int intVkey(string name, bool append=0);	// return column position if integer property with this name exists, -1 otherwise; append prototype if append==true
    int doubleVkey(string name, bool append=0);
    int stringVkey(string name, bool append=0);
    
    int intColumns() const { return intPropertyKey.size(); }
    int doubleColumns() const { return doublePropertyKey.size(); }
    int stringColumns() const { return stringPropertyKey.size(); }
    
    void setSortKey(string columnname)
    {
      intSortVkey = intVkey(columnname);
      doubleSortVkey = doubleVkey(columnname);
      if(intSortVkey<0 && doubleSortVkey<0)
	cerr << "ERROR: cannot set sort key " << columnname << endl;
    }


  static bool sortObjects(const Object* a, const Object* b)
    {
	if(a->prototype->intSortVkey >= 0)
	  return a->intProperty[a->prototype->intSortVkey] < b->intProperty[a->prototype->intSortVkey];
	if(a->prototype->doubleSortVkey >= 0)
	  return a->doubleProperty[a->prototype->doubleSortVkey] < b->doubleProperty[a->prototype->doubleSortVkey];
	cerr << "ERROR: sort key not set!" << endl;
	return 0;
    }

  private:
    int intSortVkey; int doubleSortVkey; // key of column by which this object should be sorted
    
  friend ostream& operator<<(ostream& os, ObjectPrototype& dt);
};


// class ObjectCollection:
//   a collection of objects, i.e. a catalog
//   can return filtered subsets of ObjectCollections, read stuff
//   important: if an instance of ObjectCollection reads the data itself, it will consider it 'its own' data and delete it on destruction

class ObjectCollection {
 
  public:
      
    /////////////////////////////////// CONSTRUCTORS, READING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    // create catalog
    ObjectCollection() : prototype(0), my_memory(false) { };									// creates a pretty void collection
    ObjectCollection(ObjectPrototype *dprototype) : my_memory(true) { prototype = new ObjectPrototype(dprototype); };		// use this to prepare a collection for filling in objects manually
    ObjectCollection(ifstream &file, const int keys[], const string names[], const Type types[], int n, string prefix="");	// read from ASCII catalog
      // file:   ifstream object of ASCII catalog
      // keys:   array of column numbers to be read (1,2,...)
      // names:  names of the columns for later reference
      // types:  column types (integerType, doubleType, stringType)
      // n:      number of columns to be read, i.e. the length of the previous arrays
      // prefix: prefix to prepend any of the column names given in the names vector with
    
    ObjectCollection(string   file, string extension, vector<string> columns=vector<string>(0), string prefix="");		// read from FITS_LDAC catalog
      // file:      filename
      // extension: extension name
      // columns:   vector of column names to be read; if empty, reads all
      // prefix:    prefix to prepend any of the column names with
    ObjectCollection(string   file,    int extension, vector<string> columns=vector<string>(0), string prefix="");		// read from FITS_LDAC catalog
      // file:      filename
      // extension: extension number
      // columns:   vector of column names to be read; if empty, reads all
      // prefix:    prefix to prepend any of the column names with    
    ~ObjectCollection();
    
    /////////////////////////////////// APPENDING (more lines) \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    // all these methods add new objects to the catalog
    
    void appendFile(ifstream &file); // append current catalog by reading file: i.e. more lines; this assumes the file has the same format as the one read by the constructor
    void appendFITSTable(string file, string extension, vector<string> columns); // append current catalog by reading FITS file
    
    void appendCollectionDeep(ObjectCollection &cc); 	// append (i.e. more lines) current catalog by another one; deep copy objects; this tacitly assumes the other one has the same prototype
    void makeEmpty(int nrows); 				// generate an empty catalog of nrows rows with all-zero entries
    void appendObject(Object *cc); 			// append current catalog by object; does not copy, just saves pointer

    /////////////////////////////////// EXTENDING (more columns) \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    void extendCollection(ObjectCollection &cc, string prefix=""); // extend (i.e. more columns) with a different catalog of the same objects in the same order
    void dualMatchFile(ifstream &file, string prefix, const int keys[], const string names[], const Type types[], int n, 
		       string xname, string yname, string xnamenew, string ynamenew, double r1, double r2, bool append_unmatched=true, bool shout=false); 
       // add information from additional ASCII catalog to columns with prefix "prefix"
       //  match to current objects if position xname,yname agrees within less than r1 and within less than r2 uniquely
       //  otherwise append catalog (unless you say append_unmatched, then I will forget about new lines that do not match old objects)
       //  
    void dualMatchCollection(ObjectCollection &other, string prefix, string xname, string yname, string xnamenew, string ynamenew, double r1, double r2, bool append_unmatched, bool shout=false);
       // do the same, but with a collection that maybe is already manipulated somehow
    
    /////////////////////////////////// FILTERING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    ObjectCollection* filter(Filter f) const; 		// return filtered ObjectCollection
      // note that internally no memory is copied, but rather a new vector of pointers is used to collect matching Objects; you can safely delete the new ObjectCollection or keep on filtering
    
    ObjectCollection* filter(vector<Filter> f) const; 	// return filtered ObjectCollection where filters are combined with AND  
      
    vector<unsigned char> filterFlagVector(vector<Filter> f) const; // returns vector which is 1 for good objects, 0 for deselected objects
   
    ObjectCollection* filterString(string filterkey, string filtervalue) const;  // filter for a string property
    ObjectCollection* filterString(int filterkey, string filtervalue) const;
    ObjectCollection* filterStringRemove(string filterkey, vector<string> baddata) const; // remove all objects where string property filterkey contains any of the baddata strings
    ObjectCollection* filterStringRemove(int filterkey, vector<string> baddata) const; 

    //////////////////////////////////// BOOTSTRAPPING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
    ObjectCollection* bootstrapWeight(string weightkey, double weightsum) const;
    ObjectCollection* bootstrapWeight(int weightkey, double weightsum) const;
    // return bootstrapped sample with weightkey sum equal (or slightly larger) than weightsum

    /////////////////////////////////// COLUMN OPERATIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    // simple column operations
    double sum(string colname) const;					// sum of column in all objects in the collection
    double sumProduct(string colname1, string colname2) const;		// sum of product of two columns
    double sumSq(string colname) const;					// sum of squares of column
    double sumErrSq(string col1, string col2) const; 			// sum of squared deviations between two columns -- good for chi^2
    double sumErrSqWt(string col1, string col2, string col3) const; 	// sum of squared deviations between two columns, weighted with a third one -- good for chi^2
    double sumErrSqWt(string col1, double col2, string col3) const; 	// sum of squared deviations column and model, weighted with a third one -- good for chi^2
    double sumErrSqIntObs(string col1, string col2, 			// sum sqared deviations and normalize by interr^2+obserr^2 -- good for chi^2
			  double interr, string obserr) const; 		//  note that interr is assumed constant and obserr is taken from a third column

    double average(string colname) const;				// average of a column over all objects in the collection
    double avgstdv(string colname) const;				// standard deviation of the average calculated above (goes like 1/sqrt(n-1))
    double averageWeighted(string colname, string weightname) const;    // weighted average according to weight column
    double median(string colname);					// median of column
    double percentile(string colname, double p);			// percentile of column, i.e. lowest value that is above p*100% of samples
    
    // create new or modify old columns using an arbitrary function of any number of columns
    // transformColumn changes and existing column, which is passed to the function as a first argument (there can be up to 3 additional arguments)
    // transformColumnNew creates a new column that is a function of up to 12 columns
    
    void transformColumn(string colname, double (*f)(double)); // transform column by a function; must be a double column!
    void transformColumn(string colname, double (*f)(double,double), string arg1); // transform column by a function; must be a double column!
    void transformColumn(string colname, double (*f)(double,double,double), string arg1, string arg2); // transform column by a function; must be a double column!
    void transformColumn(string colname, double (*f)(double,double,double,double), string arg1, string arg2, string arg3); // transform column by a function; must be a double column!
    int  transformColumnNew(string newcolname, double (*f)(double), string arg); // create new column from an existing one with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double), string arg1, string arg2); // create new column from two existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double), string arg1, string arg2, string arg3); // create new column from three existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double), string arg1, string arg2, string arg3, string arg4); // create new column from four existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5); // create new column from five existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6); // create new column from six existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7); // create new column from seven existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8); // create new column from eight existing ones with a function
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13, string arg14); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13, string arg14, string arg15); 
    int  transformColumnNew(string newcolname, double (*f)(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double), string arg1, string arg2, string arg3, string arg4, string arg5, string arg6, string arg7, string arg8, string arg9, string arg10, string arg11, string arg12, string arg13, string arg14, string arg15, string arg16); 

    double epanechnikovDensity(string cx, string cy, double x, double y, double r) const; // 2d epanechnikov kernel density with radius r
    
    //
    // arkane stuff that nobody will want to use...
    //
    //
    // ellipsoid statistics, something I need for beta(mag)
    double averageEllipsoid(string colname, vector<string> ecolnames,   // average of property among objects in an ellipsoid (with lost of complicated stuff I need for beta(mag)
	vector<double> centers, vector<double> radii, 
	const vector<double> maxima = vector<double>(), 
		bool talk=false) const;
    Object* randomFromEllipsoid(vector<string> ecolnames, 
		vector<double> centers, vector<double> radii) const;
    int    countInEllipsoid(vector<string> ecolnames, vector<double> centers, vector<double> radii) const;

    
    int createDoublePropertyIfNecessary(string colname, double init=-1., bool initevenifexists=true); // return column key of present or newly generated double property
    int createIntPropertyIfNecessary(string colname, int init=-1, bool initevenifexists=true);    // return column key of present or newly generated double property
    int createStringPropertyIfNecessary(string colname); // return column key of present or newly generated double property
    
    bool calculateRp(string xcol, string ycol, double x0, double y0, string rpcol="rp", string phipcol="phip", bool clip_margin=false, double z0=0.5e6); 
    bool calculateRpRADec(string racol, string deccol, double ra0, double dec0, string rpcol="rp", string phipcol="phip", double kpcperdeg=1.); 
    // calculate projected radius based on two position columns and one reference position and save in new column r and phi
    // RADec version: projected distance in angles taking into account spherical geometry
    
    bool calculateGtGx(string phicol, string g1col, string g2col, string gtcol, string gxcol);
   
    void setIntProperty(string colname, int value, int first=0, int last=-1); // set property
    void setDoubleProperty(string colname, double value, int first=0, int last=-1); // set property
    
    ///////////////////////////////////////// BINNING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    void sort(string sortcolumn); // sort objects by column
    
    ObjectCollection* bin(vector<string> averaged_columns, string binning_column, int nbins); 
      // average averaged_columns for nbins bins in binning_column, of equal population 
      // warning: may change order of this ObjectCollection; the returned table is allocated and must be disallocated by the user
    
    ObjectCollection* bin(vector<string> averaged_columns, string binning_column, vector<double> binning_limits); 
      // average averaged_columns for bins in binning_column, ranging from binning_limits[0]<=binning_column<binning_limits[1] ...
      // warning: may change order of this ObjectCollection; the returned table is allocated and must be disallocated by the user
    
    ObjectCollection* bin(vector<string> averaged_columns, vector<string> binning_columns, vector<int> nbins); 
      // average averaged_columns for nbins[i] bins in binning_columns[i], of equal population 
      // warning: may change order of this ObjectCollection; the returned table is allocated and must be disallocated by the user 
    
    void applyBinAverage(ObjectCollection *bintable, string column, vector<string> binning_columns, string binaveragecolumn_name);
    
    /////////////////////////////////// HELPFUL SMALL STUFF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      
    int size() const; 				// return number of objects
    Object* operator[] (const int i) const; 	// return i'th object 
    void print(vector<string> column, ostream& sink=std::cout, bool header=true) const; 
						// print a subset of columns

    
    
    
    
    vector<Object*> objects;
    ObjectPrototype *prototype;
    
  private:
    
    void setPrototype(const int keys[], const string names[], const Type types[], int n); // set object prototype from key / name / type vectors 
    Object* appendEmpty(); // generate a row with all-zero entries
    void appendPrototype(ObjectPrototype *otherPrototype, string prefix, string exception="");   // appends prototype and creates columns for new fields
    void readFITSTable(int rows, vector<vector<int> > &intcolumns,  vector<vector<double> > &doublecolumns,  vector<vector<string> > &stringcolumns);
    
    ObjectCollection* filterInt(string filterkey, int minInt, int maxInt) const;
    ObjectCollection* filterInt(int filterkey, int minInt, int maxInt) const;
    
    ObjectCollection* filterDouble(string filterkey, double minDouble, double maxDouble) const;
    ObjectCollection* filterDouble(int filterkey, double minDouble, double maxDouble) const;
    
    bool my_memory;
    
    friend ostream& operator<<(ostream& os, ObjectCollection& dt);
};

///// some stuff useful for filtering and column operations

// class Region (and relatives): 
//   sometimes useful for filtering regions
class Region {
public:
  virtual double is_inside(double x1)=0;
  virtual double is_inside(double x1, double x2)=0;
};

class Rectangle : public Region {
public:
  Rectangle(double dx1, double dy1, double dx2, double dy2) : x1(dx1), x2(dx2), y1(dy1), y2(dy2) { }
  
  double is_inside(double x) { cerr << "rectangle is 2d" << endl; return 0; }
  double is_inside(double x, double y) {
   if(x1<=x && x<=x2 && y1<=y && y<=y2)
     return 1;

   return 0;
  }
  
private:
  double x1,x2,y1,y2;
};

class RegionSum : public Region {
public:
  RegionSum(Region &region) : me(region), next(0) { }
  ~RegionSum() { if(next) delete next; }
  void append(Region &appendee) {
    if(next) next->append(appendee);
    else {
	next = new RegionSum(appendee);
    }
  }
  double is_inside(double x) {
    if(next) return me.is_inside(x)+next->is_inside(x);
    return me.is_inside(x);
  }
  double is_inside(double x, double y) {
    if(next) return me.is_inside(x,y)+next->is_inside(x,y);
    return me.is_inside(x,y);
  }
  
private:
  Region &me;
  RegionSum *next;
};

class InverseRegion : public Region {
public: 
  InverseRegion(Region& dinverse) : inverse(dinverse) { }
  
  double is_inside(double x) { return 1.-inverse.is_inside(x); }
  double is_inside(double x, double y) { return 1.-inverse.is_inside(x,y); }
  
private:
  Region &inverse;
};


namespace FilterFunctions {
// use these with ObjectCollection::transformColumnNew for simple column operations
static double copy(double a)
{
 return a; 
}
static double add(double a, double b)
{
 return a+b; 
}
static double subtract(double a, double b)
{
 return a-b; 
}
static double square(double a)
{
 return a*a; 
}
static double multiply(double a, double b)
{
 return a*b; 
}
static double divide(double a, double b)
{
 return a/b; 
}
static double hypot(double a, double b)
{
 return sqrt(a*a+b*b); 
}
static double ranf()
{
 return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}
static double gaussian_random(float m, float s)	/* normal random variate generator */
{				        	/* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
static double inverseSquare(double m)
{
 return 1./m/m; 
}

static double pow2(double x)
{
 return x*x; 
}

template <typename T>
static string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

}

#endif
