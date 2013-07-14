#include "PreditorPreyInterface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"

PreditorPreyInterface::PreditorPreyInterface(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        int N, double a, double b, double s)  : 
  globalData(global_data),
  initialGuess(N),
  lambda(a),
  beta(b),
  scale(s),
  n(N),
  outputFilePtr(NULL)
{
  init();
}

PreditorPreyInterface::PreditorPreyInterface(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        int N, double a, double b, double s, ofstream& file)  : 
  globalData(global_data),
  initialGuess(N),
  lambda(a),
  beta(b),
  scale(s),
  n(N),
  outputFilePtr(&file)
{
  init();
}

const NOX::LAPACK::Vector&
PreditorPreyInterface::getInitialGuess()
{
  return initialGuess;
}

bool
PreditorPreyInterface::computeF(NOX::LAPACK::Vector& f, 
             const NOX::LAPACK::Vector &x)
{
  double h = 1.0 / double(n-1);
  double hh = h*h;

  f(0) = x(0) - alpha;
  for (int i=1; i<n-1; i++)
    f(i) = D1*(x(i-1) - 2.0*x(i) + x(i+1)) / hh
      + alpha - (beta + 1.0)*x(i) + x(i)*x(i)*x(i+n);
  f(n-1) = x(n-1) - alpha;

  f(n) = x(n) - beta/alpha;
  for (int i=n+1; i<2*n-1; i++)
    f(i) = D2*(x(i-1) - 2.0*x(i) + x(i+1)) / hh
      + beta*x(i-n) - x(i-n)*x(i-n)*x(i);
  f(2*n-1) = x(2*n-1) - beta/alpha;
  
  return true;
}

bool
PreditorPreyInterface::computeJacobian(NOX::LAPACK::Matrix<double>& J, 
             const NOX::LAPACK::Vector & x)
{
  double h = 1.0 / double(n-1);
  double hh = h*h;

  J(0,0) = 1.0;
  for (int i=1; i<n-1; i++) {
    J(i,i-1) = D1/hh;
    J(i,i) = 3.0 - 6.0*x(i) -x(i+n) -5*lambda*exp( -5.0*x(i) );
    J(i,i+1) = D1/hh;
    J(i,i+n) = -x(i);
  }
  J(n-1,n-1) = 1.0;

  J(n,n) = 1.0;
  for (int i=n+1; i<2*n-1; i++) {
    J(i,i-n) = 3.0*x(i) - 2.0*x(i-n)*x(i);
    J(i,i-1) = D2/hh;
    J(i,i) = -2.0*D2/hh - x(i-n)*x(i-n);
    J(i,i+1) = D2/hh;
  }
  J(2*n-1,2*n-1) = 1.0;

  return true;
}

void
PreditorPreyInterface::setParams(const LOCA::ParameterVector& p) {
  lambda = p.getValue("lambda");
  beta = p.getValue("beta");
  scale = p.getValue("scale");
}

void
PreditorPreyInterface::init() {

  for (int i=0; i<n; i++) 
    initialGuess(i) = 
      i*(n-1-i)*source_param(lambda, scale)/((n-1)*(n-1)) + 0.001;
}

void
PreditorPreyInterface::printSolution(const NOX::LAPACK::Vector &x,
                                    const double conParam)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() 
      << "At parameter value: " << setprecision(8) << conParam
      << "   the solution vector (norm="<<x.norm()<<") is\n";

    if (n < 8) {
      for (int i=0; i<n; i++) 
  globalData->locaUtils->out() << " " << x(i);
    }
    else {
      for (int i=0; i<6; i++) 
  globalData->locaUtils->out() << " " << x(i);
       globalData->locaUtils->out() << " ...";
      for (int i=n-2; i<n; i++)  
  globalData->locaUtils->out() << " " << x(i);
    }
    globalData->locaUtils->out() << std::endl;
  }

  if (outputFilePtr != NULL) {
    (*outputFilePtr) << conParam << " ";
    for (int i=0; i<n; i++)
      (*outputFilePtr) << x(i) << " ";
    (*outputFilePtr) << std::endl;
  }
}
