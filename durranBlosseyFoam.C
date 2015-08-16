#include "fvCFD.H"
#include "simpleMatrix.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()

/*    const scalar x_min = -150e3;
    const scalar x_max = 150e3;
    const scalar dx = (x_max - x_min) / 300;
*/

    while (runTime.loop())
    {
    }

    Field<scalar> rhs(4);
    rhs[0] = 0;
    rhs[1] = 2;
    rhs[2] = 4;
    rhs[3] = 6;

    simpleMatrix<scalar> m(4, 0, 0);
    m.source() = rhs;

    for (int i=0; i<4; i++) {
        m[i][i] = 2;
    }

    const Field<scalar>& solution = m.solve();

    Info << solution << endl;

    return 0;
}

