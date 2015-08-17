#include "fvCFD.H"
#include "simpleMatrix.H"

void initialise(Field<scalar>& phi, dimensionedScalar x_min, dimensionedScalar dx);
void write(Field<scalar>& phi, dimensionedScalar t, dimensionedScalar x_min, dimensionedScalar dx);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()

    const dimensionedScalar x_min("x_min", dimLength, -150e3);
    const dimensionedScalar x_max("x_max", dimLength, 150e3);
    const int nx = 300;
    const dimensionedScalar dx("dx", (x_max - x_min) / nx);
    const dimensionedScalar u0("u0", dimVelocity, 10);
    const dimensionedScalar c("c", u0 * dt / dx);

    Info << "# Courant no " << c.value() << endl;
    Info << "# t x theta" << endl;
   
    Field<scalar> theta_old(nx);
    Field<scalar> theta(nx);
    Field<scalar> theta_new(nx);

    dimensionedScalar t("t", dimTime, 0);
    initialise(theta, x_min, dx);
    write(theta, t, x_min, dx);

    theta_old = theta;
    scalar dt_multiplier = 0.5; // forward-in-time for the first timestep

    while (runTime.loop())
    {
        forAll(theta_new, I)
        {
            theta_new[I] = theta_old[I] - dt_multiplier*c.value()*(theta[(I+1)%nx] - theta[(I-1)%nx]);
        }

        t += dt;
        dt_multiplier = 1;

        if (static_cast<int>(t.value()) % 5000 == 0)
        {
            write(theta_new, t, x_min, dx);
        }

        theta_old = theta;
        theta = theta_new;
    }

    return 0;
}

void initialise(Field<scalar>& phi, dimensionedScalar x_min, dimensionedScalar dx)
{
    dimensionedScalar x0("x0", dimLength, -50e3);
    dimensionedScalar Ax("Ax", dimLength, 25e3);
//    scalar z0 = 9e3;
//    scalar Az = 3e3;

    dimensionedScalar x = x_min + dx/2;
    forAll(phi, I)
    {
        dimensionedScalar r = sqrt(sqr((x-x0)/Ax));
        phi[I] = (r.value() <= 1) ? sqr(Foam::cos(M_PI*r.value()/2)) : 0;
        x += dx;
    }
}

void write(Field<scalar>& phi, dimensionedScalar t, dimensionedScalar x_min, dimensionedScalar dx)
{
    dimensionedScalar x = x_min + dx/2;
    forAll(phi, I)
    {
        Info << t.value() << " ";
        Info << x.value() << " " << phi[I] << endl;
        x += dx;
    }
    Info << endl << endl;
}
