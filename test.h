#ifndef TEST_H
#define TEST_H
//#include "OP.h"
#include "Sub.h"
void test();
void test()
{
        
        Sub a(1), b(2);

        Sub c(2, a, b, 0.5);

        //Sub d(3,c,a,0.5);
        c.show();c.save();
        Sub d;d.read(2);
        cout<<"haha"<<endl;
        d.show();

        /*OP a(SpinZ);
        OP b(SpinZ);

        OP c;
        c.kron(b,a);
        c.show();*/

}

#endif // TEST_H
