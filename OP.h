#ifndef OP_H
#define OP_H

#include <string>
#include <fstream>
#include <unordered_map>
#include <map>
#include <Eigen/Eigenvalues> 
#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <sstream>
#include "Parameter.h"

struct classcom
{
        size_t operator()(const std::pair<int, int>& l) const
        {
                return (size_t) ((l.first<<20) ^l.second);
        }
};

using namespace Eigen;
using namespace std;

enum spintype
{
        SpinEye, SpinZ, SpinP, SpinM
};

std::string itos(int i);

class OP
{
private:
        unordered_map<int, int> _QDim;
        map<int, MatrixXd> _QMat;
        unordered_map<int, int> _RLQ;
public:
        static int Max;
        const unordered_map<int, int>& QDim()const{return _QDim;};
        const map<int, MatrixXd>& QMat()const{return _QMat;};
        const unordered_map<int, int>& RLQ()const{return _RLQ;};

        OP(){};
        ~OP(){};
        OP(const spintype& type);
        void Initial(const spintype& type);
        OP(const OP& a);

        void findDim(const OP& a, const OP& b, unordered_map<int, int>& oldDim, unordered_map<std::pair<int, int>, int, classcom>& startDim) const;
        void kron(const OP& a, const OP&b);
        void trans(const OP& a);

        void add(const OP& a, const OP& b);               //for the operator +, but don't need to copy anything.
        void add(const OP& a);                 //fot the operator +=, but don't need to copy anything.
        void time(const double& x);
        void time(const OP& a, const double& x);
        void time(const double& x, const OP& a);
        void time(const OP& a, const OP& b);
        OP operator+(const OP& a);
        OP operator*(const double& x);
        OP operator*(const OP& a);
        void operator=(const OP& a)
        {       
                clear();
                _QDim = a._QDim;
                _RLQ = a._RLQ;
                _QMat = a._QMat;
        };


        int ltime(const OP& a, const OP& wave);      
        //this function is used for the operator a function on the wave operator wave. 
        //the return of this function is used to see whether the time is finished or not.
        int rtime(const OP& a, const OP& wave);      
        //this two functions are used for the operator of System(the a in l), and Envirment(the a in r)

        //wave operator, for which don't have the QDim. (the operator in the QWave can't have a QDim)
        void addWave(const OP& a, const OP& b);
        void addWave(const OP& a);


        //for the truncation
        void getTruncU(const Parameter& para, const OP& OPWave);
        void getTruncUR(const Parameter& para, const OP& OPWave);
        void getTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr);
        void getTruncUR(const Parameter& para, const OP& OPWave, double& trance, double& truncerr);
        void truncL(const OP& trunc, const OP& O);
        void truncR(const OP& O, const OP& trunc);
        void trunc(const OP& truncU);


        void getDenS(const OP& OPWave);
        void getDenE(const OP& OPWave);
        void DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr);
        void DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr, double& Entanglement);


        void save(std::ofstream& outfile)const;
        void read(std::ifstream& infile);



        void truncsave(const int& orbital);
        void truncread(const int& orbital);




        void show()const;
        void clear()
        {
                _RLQ.clear();
                _QMat.clear();
                _QDim.clear();
        };

        friend class QWave;

};


#endif // OP_H
