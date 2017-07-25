#include "OP.h"

struct Eigstruct
{
        int q;
        double lamda;
        VectorXd state;
};


std::string itos(int i)
{
        std::stringstream s;
        s << i;
        return s.str();
}

bool comp(const Eigstruct& a, const Eigstruct& b);
bool comp(const Eigstruct& a, const Eigstruct& b)
{
        return (a.lamda > b.lamda);
}




OP::OP(const spintype& type)
{
        for(int i=0; i<=2; ++i)
        {
                _QDim.insert(pair<int, int>(i, 1));
        }

        switch (type)
        {
                case SpinEye:
                {
                        for(int i=0; i<=2; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i));
                                MatrixXd temp(1,1);
                                temp<<1;
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
                case SpinZ:
                {
                        for(int i=0; i<=2; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i));
                                MatrixXd temp(1,1);
                                temp<<i-1;
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
                case SpinP:
                {
                        for(int i=0; i<=1; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i+1));
                                MatrixXd temp(1,1);
                                temp<<sqrt(2);
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
                case SpinM:
                {
                        for(int i=1; i<=2; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i-1));
                                MatrixXd temp(1,1);
                                temp<<sqrt(2);
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
        }
}


void OP::Initial(const spintype& type)
{

        clear();

        for(int i=0; i<=2; ++i)
        {
                _QDim.insert(pair<int, int>(i, 1));
        }

        switch (type)
        {
                case SpinEye:
                {
                        for(int i=0; i<=2; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i));
                                MatrixXd temp(1,1);
                                temp<<1;
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
                case SpinZ:
                {
                        for(int i=0; i<=2; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i));
                                MatrixXd temp(1,1);
                                temp<<i-1;
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
                case SpinP:
                {
                        for(int i=0; i<=1; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i+1));
                                MatrixXd temp(1,1);
                                temp<<sqrt(2);
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
                case SpinM:
                {
                        for(int i=1; i<=2; ++i)
                        {
                                _RLQ.insert(pair<int, int>(i, i-1));
                                MatrixXd temp(1,1);
                                temp<<sqrt(2);
                                _QMat.insert(pair<int, MatrixXd>(i, temp));
                        }
                        break;
                }
        }
}

OP::OP(const OP& a):
_QDim(a._QDim),
_QMat(a._QMat),
_RLQ(a._RLQ)
{

}


void OP::show()const
{
        
        cout<<"==========The dimension:=========="<<endl;
        for(auto it=_QDim.begin(); it!=_QDim.end(); ++it)
        {
                cout<<it->first<<"  =>  "<<it->second<<endl;
        }
        cout<<"==========The RLQ:=========="<<endl;
        for(auto it=_RLQ.begin(); it!=_RLQ.end(); ++it)
        {
                cout<<it->first<<"  =>  "<<it->second<<endl;
        }
        cout<<"==========The Matrix:=========="<<endl;
        for(auto it=_QMat.begin(); it!=_QMat.end(); ++it)
        {
                cout<<it->first<<"  =>"<<endl<<it->second<<endl;
        }
}


void OP::findDim(const OP& a, const OP& b, std::unordered_map<int, int> &newDim, std::unordered_map<std::pair<int, int>, int, classcom> &startDim) const
{
        for (int na = 0; na <= OP::Max; ++na)
        {
                auto ita = a._QDim.find(na);
                if(ita == a._QDim.end()) continue;
                for (int nb = 0; nb <= OP::Max; ++nb)
                {
                        auto itb = b._QDim.find(nb);
                        if(itb == b._QDim.end()) continue;
                        int tempQ;
                        tempQ=ita->first + itb->first;

                        bool i = (tempQ <= Max);
                        //bool j = ((a.RLQ.at(ita->first).QN + b.RLQ.at(itb->first).QN) <= Max);

                        if (i)
                        {
                                auto  itc = newDim.find(tempQ);
                                if (itc == newDim.end())
                                {
                                        newDim.insert(std::pair<int, int>(tempQ, ita->second * itb->second));
                                        startDim.insert(pair<pair<int, int>, int>(pair<int, int>(ita->first, itb->first),0));
                                }
                                else
                                {
                                        startDim.insert(pair<pair<int, int>, int>(pair<int, int>(ita->first, itb->first),newDim[itc->first]));

                                        //startDim[std::pair<int, int>(ita->first, itb->first)] = newDim[itc->first];
                                        newDim.at(itc->first) += ita->second * itb->second;
                                }
                        }
                }

        }
}


void OP::kron(const OP &a, const OP &b)
{
        clear();





        //find the dimention of each good quantum number, and label the place to put the kron of two blocks.
        unordered_map<std::pair<int, int>, int, classcom> startDim;
        findDim(a, b, _QDim, startDim);



        //first calculate the good quantum number.
        for (auto ita = a._RLQ.begin(); ita != a._RLQ.end(); ita++)
        {
                for (auto itb = b._RLQ.begin(); itb != b._RLQ.end(); itb++)
                {
                        int tempQR;
                        tempQR = ita->first + itb->first;

                        bool i = (tempQR <= Max);
                        bool j = ((a._RLQ.at(ita->first) + b._RLQ.at(itb->first)) <= Max);

                        if (i&&j)
                        {



                                auto itc = _RLQ.find(tempQR);
                                if (itc == _RLQ.end())
                                {
                                        int tempQL;
                                        tempQL = ita->second + itb->second;
                                        _RLQ.insert(std::pair<int, int>(tempQR, tempQL));



                                        int col = _QDim.at(tempQR);
                                        int row = _QDim.at(tempQL);

                                        //creat a zero mat matched the dimention of good quantum number.
                                        MatrixXd tempm(MatrixXd::Zero(row, col));
                                        
                                        _QMat.insert(std::pair<int, MatrixXd>(tempQR, tempm));
                                }

                        }
                }
        }









        //calculate the kron and put it in the place.
        for (auto ita = a._QMat.begin(); ita != a._QMat.end(); ita++)
        {
                for (auto itb = b._QMat.begin(); itb != b._QMat.end(); itb++)
                {




                        int tempQR;
                        tempQR = ita->first + itb->first;

                        bool i = (tempQR <= Max);
                        bool j = ((a._RLQ.at(ita->first) + b._RLQ.at(itb->first) )<= Max);

                        if (i&&j)
                        {


                                

                                


                                int startR = startDim[std::pair<int, int>(ita->first, itb->first)];




                                int startL = startDim[std::pair<int, int>(a._RLQ.at(ita->first), b._RLQ.at(itb->first))];



                                for (int m = 0; m < ita->second.rows(); m++)
                                {
                                        for (int n = 0; n < ita->second.cols(); n++)
                                        {
                                                int a = itb->second.rows();
                                                int b = itb->second.cols();
                                                
                                                //QMat.at(tempQR).block(startL + m*itb->second.rows(), startR + n*itb->second.cols(),a,b)
                                                //      = ita->second(m, n)*itb->second;
                                                for (int i = 0; i < itb->second.rows(); ++i)
                                                {
                                                        for (int j = 0; j < itb->second.cols(); ++j)
                                                        {
                                                                _QMat.at(tempQR)(startL + m*itb->second.rows() + i, startR + n*itb->second.cols() + j)
                                                                                = ita->second(m, n)*itb->second(i,j);
                                                        }
                                                }
                                                
                                        }
                                }



                                
                        }

                }
        }
}


void OP::trans(const OP& a)
{
        _QDim = a._QDim;

        for (auto it = a._RLQ.begin(); it != a._RLQ.end(); it++)
        {
                _RLQ.insert(std::pair<int,int>(it->second, it->first));
                _QMat.insert(std::pair<int, MatrixXd>(it->second, a._QMat.at(it->first).transpose()));
        }

}



void OP::add(const OP&a, const OP& b)
{

        if(a._QDim.size()==0)
        {
                _QDim=b._QDim;
                _QMat=b._QMat;
                _RLQ=b._RLQ;
                return;
        }else if(b._QDim.size()==0)
        {
                _QDim=a._QDim;
                _QMat=a._QMat;
                _RLQ=a._RLQ;
                return;
        }

        if (a._QDim.size() != b._QDim.size())
        {
                cerr << "the Q dim size doesn't match!" << endl;
                exit(1);
        }

        for (auto tempit = a._QDim.begin(); tempit != a._QDim.end(); tempit++)
        {
                if (tempit->second != b._QDim.at(tempit->first))
                {
                        cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }
        }

        _QDim = a._QDim;
        _RLQ = a._RLQ;
        _QMat = a._QMat;





        for (auto tempit = b._QMat.begin(); tempit != b._QMat.end(); tempit++)
        {
                auto it = _QMat.find(tempit->first);
                if (it == _QMat.end())
                {
                        _QMat.insert(std::pair<int, MatrixXd>(tempit->first, b._QMat.at(tempit->first)));
                        _RLQ.insert(std::pair<int, int>(tempit->first, b._RLQ.at(tempit->first)));

                }
                else
                {
                        if (a._RLQ.at(tempit->first) != b._RLQ.at(tempit->first))
                        {
                                cerr << "add0: RQ != LQ, WRONG!" << std::endl;
                                exit(1);
                        }
                        _QMat.at(tempit->first) += b._QMat.at(tempit->first);
                }
        }



}


void OP::add(const OP& a)
{
        if(_QDim.size()==0)
        {
                _QDim=a._QDim;
                _QMat=a._QMat;
                _RLQ=a._RLQ;
                return;
        }else if(a._QDim.size()==0)
        {
                return;
        }

        if (_QDim.size() != a._QDim.size())
        {
                cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }

        for (auto tempit = _QDim.begin(); tempit != _QDim.end(); tempit++)
        {
                if (tempit->second != a._QDim.at(tempit->first))
                {
                        cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }
        }
        for (auto tempit = a._QMat.begin(); tempit != a._QMat.end(); tempit++)
        {
                auto it = _QMat.find(tempit->first);
                if (it == _QMat.end())
                {
                        _QMat.insert(std::pair<int, MatrixXd>(tempit->first, a._QMat.at(tempit->first)));
                        _RLQ.insert(std::pair<int, int>(tempit->first, a._RLQ.at(tempit->first)));
                }

                else
                {
                        /*if(RLQ.at(tempit->first).QN != a.RLQ.at(tempit->first).QN)
                        {
                        std::cout<<"add1: RQ != LQ, WRONG!"<<std::endl;
                        exit(1);
                        }*/

                        _QMat.at(tempit->first) += a._QMat.at(tempit->first);

                }
        }
}



OP OP::operator +(const OP& a)
{
        OP sum;


        //test the QDim.
        if (_QDim.size() != a._QDim.size())
        {
                cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }
        for (auto tempit = _QDim.begin(); tempit != _QDim.end(); tempit++)
        {
                if (tempit->second != a._QDim.at(tempit->first))
                {
                        cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }

        }



        //test the RLQ.
        

        sum._QDim = _QDim;
        sum._RLQ = _RLQ;
        sum._QMat = _QMat;




        for (auto tempit = a._QMat.begin(); tempit != a._QMat.end(); tempit++)
        {
                auto it = _QMat.find(tempit->first);
                if (it == _QMat.end())
                {
                        sum._QMat.insert(std::pair<int, MatrixXd>(tempit->first, a._QMat.at(tempit->first)));
                        sum._RLQ.insert(std::pair<int, int>(tempit->first, a._RLQ.at(tempit->first)));
                }

                else
                {
                        if (_RLQ.at(tempit->first) != a._RLQ.at(tempit->first))
                        {
                                std::cout << "+: RQ != LQ, WRONG!" << std::endl;
                                exit(1);
                        }
                        sum._QMat.at(tempit->first) += a._QMat.at(tempit->first);
                }
        }




        return sum;
}

void OP::time(const double& x)
{
        for (auto it = _QMat.begin(); it != _QMat.end(); it++)
        {
                _QMat.at(it->first) = _QMat.at(it->first) * x;
        }
}


void OP::time(const OP& a, const double& x)
{
        clear();
        _QDim = a._QDim;
        _RLQ = a._RLQ;
        for (auto it = a._QMat.begin(); it != a._QMat.end(); it++)
        {
                _QMat[it->first] = a._QMat.at(it->first)*x;
        }
}


void OP::time(const double& x, const OP& a)
{
        clear();
        _QDim = a._QDim;
        _RLQ = a._RLQ;
        for (auto it = a._QMat.begin(); it != a._QMat.end(); it++)
        {
                _QMat[it->first] = a._QMat.at(it->first)*x;
        }
}

OP OP::operator*(const double& x)
{
        OP times;


        times._QDim = _QDim;
        times._RLQ = _RLQ;


        //the matrix times the parameter number.
        for (auto tempit = _QMat.begin(); tempit != _QMat.end(); tempit++)
        {
                times._QMat.insert(std::pair<int, MatrixXd>(tempit->first, x * _QMat.at(tempit->first)));
        }

        return times;
}

OP operator*(const double& x, OP& a)
{
        return a*x;
}



void OP::time(const OP&a, const OP&b)
{
        if (a._QDim.size() != b._QDim.size())
        {
                cerr << "the Q dim size doesn't match!" << std::endl;
                exit(1);
        }
        for (auto tempit = a._QDim.begin(); tempit != a._QDim.end(); tempit++)
        {
                if (tempit->second != b._QDim.at(tempit->first))
                {
                        cerr << "the Q dimention doesn't match each other!" << std::endl;
                        exit(1);
                }

        }

        _QDim = a._QDim;




        //for the times of two operator, it is a good block systerm. just cross it.
        for (auto ita = b._QMat.begin(); ita != b._QMat.end(); ita++)
        {
                auto it = a._QMat.find(b._RLQ.at(ita->first));
                if (it != a._QMat.end())
                {

                        _RLQ.insert(std::pair<int, int>(ita->first, a._RLQ.at(it->first)));
                        _QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

                }
        }

}



int OP::ltime(const OP& a, const OP& wave)
{
        clear();
        int flag(0);
        for (auto ita = wave._QMat.begin(); ita != wave._QMat.end(); ita++)
        {
                auto it = a._QMat.find(wave._RLQ.at(ita->first));
                if (it != a._QMat.end())
                {
                        flag++;
                        _RLQ.insert(std::pair<int, int>(ita->first, a._RLQ.at(it->first)));
                        _QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second * ita->second));

                }
        }
        return flag;

}



int OP::rtime(const OP& a, const OP& wave)
{
        clear();
        int flag(0);

        for (auto ita = a._QMat.begin(); ita != a._QMat.end(); ita++)
        {
                for (auto itQ = wave._RLQ.begin(); itQ != wave._RLQ.end(); itQ++)
                {

                        if (ita->first == itQ->first)
                        {
                                flag++;
                                int tempQ(a._RLQ.at(ita->first));
                                _RLQ[tempQ] = itQ->second;
                                _QMat[tempQ] = (wave._QMat.at(itQ->first)) * ita->second.transpose();

                        }
                }
        }
        return flag;
}


void OP::addWave(const OP& a, const OP& b)
{
        clear();
        _QMat = a._QMat;
        _RLQ = a._RLQ;
        for (auto it = b._RLQ.begin(); it != b._RLQ.end(); it++)
        {
                auto tempit = a._RLQ.find(it->first);
                if (tempit != a._RLQ.end())
                {
                        if (it->second == tempit->second)
                        {
                                _QMat.at(tempit->first) += b._QMat.at(it->first);
                        }
                        else
                        {
                                cerr << "addWave: a.RLQ.second != b.RLQ.second" << std::endl;
                                exit(1);
                        }
                }
                else
                {
                        _RLQ.insert(std::pair<int, int>(it->first, it->second));
                        _QMat.insert(std::pair<int, MatrixXd>(it->first, b._QMat.at(it->first)));
                }
        }
}


void OP::addWave(const OP& a)
{
        for (auto it = a._RLQ.begin(); it != a._RLQ.end(); it++)
        {
                auto tempit = _RLQ.find(it->first);

                if (tempit != _RLQ.end())
                {
                        _QMat.at(tempit->first) += a._QMat.at(it->first);
                }
                else
                {
                        _RLQ.insert(std::pair<int, int>(it->first, it->second));
                        _QMat.insert(std::pair<int, MatrixXd>(it->first, a._QMat.at(it->first)));
                }
        }
}


//================for the truncate=======================
void OP::getTruncU(const Parameter& para, const OP& OPWave)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
                for (int i = 0; i<svd.singularValues().size(); i++)
                {
                        Eigstruct temp = { OPWave._RLQ.at(it->first), svd.singularValues()(i), svd.matrixU().col(i) };
                        Denmat.push_back(temp);
                }
        }
        //std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
        sort(Denmat.begin(), Denmat.end(), comp);


        int min = Denmat.size() < para.D() ? Denmat.size() : para.D();
        //std::cout<<"min = "<<min<<std::endl;
        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = _QDim.find(Denmat.at(i).q);
                if (itt != _QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        _QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        _RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                _QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                _QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }
        
}



void OP::getTruncUR(const Parameter& para, const OP& OPWave)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
                for (int i = 0; i<svd.singularValues().size(); i++)
                {
                        Eigstruct temp = { it->first, svd.singularValues()(i), svd.matrixV().col(i) };
                        Denmat.push_back(temp);
                }
        }
        //std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
        sort(Denmat.begin(), Denmat.end(), comp);
        //if(Denmat.size() < para.D) return 0;
        int min = (Denmat.size() < para.D()) ? Denmat.size() : para.D();
        //std::cout<<"min = "<<min<<std::endl;
        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = _QDim.find(Denmat.at(i).q);
                if (itt != _QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        _QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        _RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                _QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                _QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }

        //return 1;
}


void OP::getTruncUR(const Parameter& para, const OP& OPWave, double& trance, double& truncerr)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
                for (int i = 0; i<svd.singularValues().size(); i++)
                {
                        Eigstruct temp = { it->first, svd.singularValues()(i), svd.matrixV().col(i) };
                        Denmat.push_back(temp);
                }
        }
        //std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
        sort(Denmat.begin(), Denmat.end(), comp);
        //if(Denmat.size() < para.D) return 0;
        int min = (Denmat.size() < para.D()) ? Denmat.size() : para.D();
        //std::cout<<"min = "<<min<<std::endl;
        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = _QDim.find(Denmat.at(i).q);
                if (itt != _QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        _QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        _RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                _QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                _QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }

        //return 1;


        trance = 0;
        truncerr = 0;
        for (int i = 0; i<Denmat.size(); i++)
        {
                trance += (Denmat.at(i).lamda)*(Denmat.at(i).lamda);
        }
        for (int i = 0; i< min; i++)
        {
                truncerr += Denmat.at(i).lamda*(Denmat.at(i).lamda);
        }

        truncerr = trance - truncerr;
}


void OP::getTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                JacobiSVD<MatrixXd> svd(it->second, ComputeThinU | ComputeThinV);
                for (int i = 0; i<svd.singularValues().size(); i++)
                {
                        Eigstruct temp = { OPWave._RLQ.at(it->first), svd.singularValues()(i), svd.matrixU().col(i) };
                        Denmat.push_back(temp);
                }
        }
        //std::cout<<"Denmat.size="<<Denmat.size()<<std::endl;
        sort(Denmat.begin(), Denmat.end(), comp);


        int min = Denmat.size() < para.D() ? Denmat.size() : para.D();
        //std::cout<<"min = "<<min<<std::endl;
        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = _QDim.find(Denmat.at(i).q);
                if (itt != _QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        _QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        _RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                _QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                _QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }


        trance = 0;
        truncerr = 0;
        for (int i = 0; i<Denmat.size(); i++)
        {
                trance += pow(Denmat.at(i).lamda, 2);
        }
        for (int i = 0; i< min; i++)
        {
                truncerr += pow(Denmat.at(i).lamda, 2);
        }

        truncerr = trance - truncerr;

}



void OP::DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave.QMat().begin(); it != OPWave.QMat().end(); it++)
        {
                SelfAdjointEigenSolver<MatrixXd> es(it->second);
                for (int i = 0; i<es.eigenvalues().size(); i++)
                {
                        Eigstruct temp = { it->first, es.eigenvalues()(i), es.eigenvectors().col(i) };
                        Denmat.push_back(temp);
                }
        }

        std::stable_sort(Denmat.begin(), Denmat.end(), comp);

        int min = (Denmat.size() < para.D()) ? Denmat.size() : para.D();

        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = _QDim.find(Denmat.at(i).q);
                if (itt != _QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        _QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        _RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                _QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                _QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }
        trance = 0;
        truncerr = 0;
        for (int i = 0; i<Denmat.size(); i++)
        {
                trance += Denmat.at(i).lamda;
        }
        for (int i = 0; i< min; i++)
        {
                truncerr += Denmat.at(i).lamda;
        }

        truncerr = trance - truncerr;
}





void OP::DengetTruncU(const Parameter& para, const OP& OPWave, double& trance, double& truncerr, double& Entanglement)
{
        clear();
        std::vector<Eigstruct> Denmat;

        //get the denmat
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                SelfAdjointEigenSolver<MatrixXd> es(it->second);
                for (int i = 0; i<es.eigenvalues().size(); i++)
                {
                        Eigstruct temp = { it->first, es.eigenvalues()(i), es.eigenvectors().col(i) };
                        Denmat.push_back(temp);
                }
        }

        std::stable_sort(Denmat.begin(), Denmat.end(), comp);

        int min = (Denmat.size() < para.D()) ? Denmat.size() : para.D();

        //get the RLQ/QDim
        for (int i = 0; i<min; i++)
        {
                auto itt = _QDim.find(Denmat.at(i).q);
                if (itt != _QDim.end())
                {
                        itt->second += 1;
                }
                else
                {
                        _QDim.insert(std::pair<int, int>(Denmat.at(i).q, 1));
                        _RLQ.insert(std::pair<int, int>(Denmat.at(i).q, Denmat.at(i).q));
                }
        }

        //get the QMat
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                int L = Denmat.at(i).state.size();
                                int R = it->second;
                                MatrixXd tempmat(L, R);
                                _QMat.insert(std::pair<int, MatrixXd>(it->first, tempmat));
                                break;
                        }
                }
        }

        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                int ord(0);
                for (int i = 0; i<min; i++)
                {
                        if (Denmat.at(i).q == it->first)
                        {
                                _QMat.at(it->first).col(ord) = Denmat.at(i).state;
                                ord++;
                        }
                }
        }
        trance = 0;
        truncerr = 0;
        Entanglement = 0;
        for (int i = 0; i<Denmat.size(); i++)
        {
                trance += Denmat.at(i).lamda;
                if(Denmat.at(i).lamda > 0.0000000001)
                Entanglement = Entanglement - Denmat.at(i).lamda*log(Denmat.at(i).lamda)/log(2);
                
        }
        for (int i = 0; i< min; i++)
        {
                truncerr += Denmat.at(i).lamda;
        }

        truncerr = trance - truncerr;
}



void OP::truncL(const OP& trunc, const OP& O)
{
        clear();
        if (O._QDim.size() == 0) exit(1);

        for (auto ita = O._QDim.begin(); ita != O._QDim.end(); ita++)
        {
                auto it = trunc._QDim.find(ita->first);
                if (it != trunc._QDim.end())
                {
                        _QDim.insert(std::pair<int, int>(ita->first, trunc._QDim.at(ita->first)));
                }
        }


        for (auto ita = O._QMat.begin(); ita != O._QMat.end(); ita++)
        {
                auto it = trunc._QMat.find(O._RLQ.at(ita->first));
                if (it != trunc._QMat.end())
                {

                        _RLQ.insert(std::pair<int, int>(ita->first, it->first));
                        _QMat.insert(std::pair<int, MatrixXd>(ita->first, it->second.transpose() * ita->second));
                }
        }
}



void OP::truncR(const OP& O, const OP& trunc)
{
        clear();

        for (auto ita = O._QDim.begin(); ita != O._QDim.end(); ita++)
        {

                auto it = trunc._QDim.find(ita->first);

                if (it != trunc._QDim.end())
                {

                        _QDim.insert(std::pair<int, int>(ita->first, trunc._QDim.at(ita->first)));

                }
        }

        for (auto ita = O._QMat.begin(); ita != O._QMat.end(); ita++)
        {

                auto it = trunc._QMat.find(ita->first);

                if (it != trunc._QMat.end())
                {

                        _RLQ.insert(std::pair<int, int>(ita->first, O._RLQ.at(ita->first)));
                        
                        _QMat.insert(std::pair<int, MatrixXd>(ita->first, ita->second * it->second));
                }
        }
}


void OP::trunc(const OP& truncU)
{
        OP temp;

        temp.truncR(*this, truncU);

        truncL(truncU, temp);
}


void OP::getDenS(const OP& OPWave)
{
        clear();
        for (auto it = OPWave._RLQ.begin(); it != OPWave._RLQ.end(); it++)
                _RLQ.insert(std::pair<int, int>(it->second, it->second));
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                MatrixXd temp = it->second * (it->second).transpose();
                _QMat.insert(std::pair<int, MatrixXd>(OPWave._RLQ.at(it->first), temp));
        }
}



void OP::getDenE(const OP& OPWave)
{
        clear();
        for (auto it = OPWave._RLQ.begin(); it != OPWave._RLQ.end(); it++)
                _RLQ.insert(std::pair<int, int>(it->first, it->first));
        for (auto it = OPWave._QMat.begin(); it != OPWave._QMat.end(); it++)
        {
                MatrixXd temp = (it->second).transpose() * it->second;
                _QMat.insert(std::pair<int, MatrixXd>(it->first, temp));
        }
}



//save the operator data.
void OP::save(std::ofstream& outfile)const
{
        //std::ofstream outfile(str);
        //save the QDim.

        outfile.precision(30);


        outfile << _QDim.size() << std::endl;
        for (auto it = _QDim.begin(); it != _QDim.end(); it++)
        {
                outfile << it->first << "       " << it->second << "        ";
        }


        //save the RLQ.
        outfile << _RLQ.size() << std::endl;
        for (auto it = _RLQ.begin(); it != _RLQ.end(); it++)
        {
                outfile << it->first << "          " << it->second << "             ";
        }
        outfile << std::endl;



        //save the QMat.
        for (auto it = _QMat.begin(); it != _QMat.end(); it++)
        {
                outfile << it->first << std::endl;

                outfile<< it->second.rows() << "     " <<it->second.cols()<<std::endl;

                outfile << it->second << std::endl;
                /*for(int i=0; i<it->second.n_rows; i++)
                {
                for(int j=0; j<it->second.n_cols; j++)
                {
                outfile <<std::setprecision(20)<< it->second(i,j)<<std::endl;
                }
                }*/

        }


        //outfile.close();
}


//read operator data from the file.
void OP::read(std::ifstream& infile)
{
        //std::ifstream infile(str);
        /*if(!infile.is_open())
        {
        std::cout<<"the file "<< str << " is not exit " << std::endl;
        exit(1);
        }*/

        //read in the QDim.
        clear();
        int size1;
        infile >> size1;

        int tempQ;
        int tempint;
        for (int i = 0; i < size1; i++)
        {

                infile >> tempQ >> tempint;
                _QDim[tempQ] = tempint;
        }

        //read in the RLQ.
        int size2;
        infile >> size2;
        int tempQ1;
        for (int i = 0; i < size2; i++)
        {
                infile >> tempQ >> tempQ1;
                _RLQ[tempQ] = tempQ1;
        }




        //read in the QMat.
        for (int it = 0; it < _RLQ.size(); ++it)
        {

                int tempQR;

                infile >> tempQR;

                int row, col;

                infile >> row >> col;

                MatrixXd A(row, col);
                for (int i = 0; i < row; i++)
                {
                        for (int j = 0; j < col; j++)
                        {
                                infile >> A(i, j);

                        }
                }

                _QMat[tempQR] = A;
        }
        //infile.close();



}



void OP::truncsave(const int& orbital)
{
        //std::ofstream outfile(str);
        //save the QDim.
        std::string str = itos(orbital);

        str = "./data/trunc" +str;


        std::ofstream outfile(str);

        outfile.precision(30);


        outfile << QDim().size() << std::endl;
        for (auto it = QDim().begin(); it != QDim().end(); it++)
        {
                outfile << it->first << "       " << it->second << "        "<<std::endl;
        }


        //save the RLQ.
        outfile << RLQ().size() << std::endl;
        for (auto it = RLQ().begin(); it != RLQ().end(); it++)
        {
                outfile << it->first << "          " << it->second << "             "<<std::endl;
        }
        outfile << std::endl;



        //save the QMat.
        for (auto it = QMat().begin(); it != QMat().end(); it++)
        {
                outfile << it->first << std::endl;

                outfile << it->second.rows()<<std::endl;

                outfile << it->second.cols() << std::endl;

                outfile << it->second << std::endl;
                

        }


        outfile.close();
}




//read the truncation operator data from the file.
void OP::truncread(const int& orbital)
{
        std::string str = itos(orbital);

        str = "./data/trunc" +str;

        std::ifstream infile(str);
        int size1;
        infile >> size1;

        int tempQ;
        int tempint;
        for (int i = 0; i < size1; i++)
        {

                infile >> tempQ >> tempint;
                _QDim[tempQ] = tempint;
        }

        //read in the RLQ.
        int size2;
        infile >> size2;
        int tempQ1;
        for (int i = 0; i < size2; i++)
        {
                infile >> tempQ >> tempQ1;
                _RLQ[tempQ] = tempQ1;
        }




        //read in the QMat.
        for (int it = 0; it < _RLQ.size(); ++it)
        {

                int tempQR;

                infile >> tempQR;
                int dimL, dimR;
                infile >>dimL >>dimR;

                MatrixXd A(dimL, dimR);
                for (int i = 0; i < dimL; i++)
                {
                        for (int j = 0; j < dimR; j++)
                        {
                                infile >> A(i, j);

                        }
                }

                _QMat[tempQR] = A;
        }
        infile.close();



}