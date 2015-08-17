#include "fvCFD.H"
#include "simpleMatrix.H"

typedef struct
{
    const dimensionedScalar& x_min;
    const dimensionedScalar& x_max;
    const dimensionedScalar& z_min;
    const dimensionedScalar& z_max;
    const int nx;
    const int nz;
    const dimensionedScalar& dx;
    const dimensionedScalar& dz;
} geometry;

void initialise(FieldField<Field, scalar>& phi, const geometry& mesh);
void write(const FieldField<Field, scalar>& phi, dimensionedScalar t, const geometry& mesh);
scalar grad_x(const FieldField<Field, scalar>& phi, const label I, const label K, const geometry& mesh);
scalar grad_z(const FieldField<Field, scalar>& phi, const label I, const label K, const geometry& mesh);
int modulo(int a, int b);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #define dt runTime.deltaT()

    geometry mesh = {
       dimensionedScalar("x_min", dimLength, -150e3),
       dimensionedScalar("x_min", dimLength, 150e3),
       dimensionedScalar("z_min", dimLength, 0),
       dimensionedScalar("z_max", dimLength, 25e3),
       300,
       50,
       dimensionedScalar("dx", (mesh.x_max - mesh.x_min)/mesh.nx),
       dimensionedScalar("dz", (mesh.z_max - mesh.z_min)/mesh.nz),
    };

    const dimensionedVector u0("u0", dimVelocity, vector(10,0,0.5));
    const scalar cx = (u0[0] * dt / mesh.dx).value();
    const scalar cz = (u0[2] * dt / mesh.dz).value();

    Info << "# t x z theta" << endl;
   
    FieldField<Field,scalar> theta(mesh.nz+1);
    FieldField<Field,scalar> theta_new(mesh.nz+1);

    dimensionedScalar t("t", dimTime, 0);
    initialise(theta, mesh);
    initialise(theta_new, mesh);
    write(theta, t, mesh);

    while (runTime.loop())
    {
        Info << "# " << runTime.timeName() << endl;
        for (int corr=0; corr < 3; corr++)
        {
            forAll(theta_new, K)
            {
                forAll(theta_new[K], I)
                {
                    theta_new[K][I] = theta[K][I] - 0.5*(
                        cx*(grad_x(theta_new, I, K, mesh) + grad_x(theta, I, K, mesh)) + 
                        cz*(grad_z(theta_new, I, K, mesh) + grad_z(theta, I, K, mesh)));
                }
            }
        }

        t += dt;

        if (static_cast<int>(t.value()) % 5000 == 0)
        {
            write(theta_new, t, mesh);
        }

        theta = theta_new;
    }

    return 0;
}

scalar grad_x(const FieldField<Field, scalar>& phi, const label I, const label K, const geometry& mesh)
{
    return 0.5*(phi[K][modulo(I+1, mesh.nx)] - phi[K][modulo(I-1, mesh.nx)]);
}

scalar grad_z(const FieldField<Field, scalar>& phi, const label I, const label K, const geometry& mesh)
{
    return 0.5*(phi[modulo(K+1, mesh.nz+1)][I] - phi[modulo(K-1, mesh.nz+1)][I]);
}

int modulo(int a, int b)
{
    return ((a % b) + b) % b;
}

void initialise(FieldField<Field,scalar>& phi, const geometry& mesh)
{
    dimensionedScalar x0("x0", dimLength, -50e3);
    dimensionedScalar z0("z0", dimLength, 9e3);
    dimensionedScalar Ax("Ax", dimLength, 25e3);
    dimensionedScalar Az("Az", dimLength, 3e3);

    dimensionedScalar z = mesh.z_min;
    forAll(phi, K)
    {
        dimensionedScalar x = mesh.x_min + mesh.dx/2;
        phi.set(K, new Field<scalar>(mesh.nx));
        forAll(phi[K], I)
        {
            dimensionedScalar r = sqrt(sqr((x-x0)/Ax) + sqr((z-z0)/Az));
            phi[K][I] = (r.value() <= 1) ? sqr(Foam::cos(M_PI*r.value()/2)) : 0;
            x += mesh.dx;
        }
        z += mesh.dz;
    }
}

void write(const FieldField<Field,scalar>& phi, dimensionedScalar t, const geometry& mesh)
{
    dimensionedScalar z = mesh.z_min;
    forAll(phi, K)
    {
        dimensionedScalar x = mesh.x_min + mesh.dx/2;
        forAll(phi[K], I)
        {
            Info << t.value() << " ";
            Info << x.value() << " " << z.value() << " " << phi[K][I] << endl;
            x += mesh.dx;
        }
        z += mesh.dz;
    }
    Info << endl << endl;
}
