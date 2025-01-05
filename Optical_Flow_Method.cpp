#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath> 
#include <sstream>
#include <chrono>
#include <string>

using namespace std;
using namespace chrono;

const int sizes(100);
const double finalTime(0.05);
const double timeStep(0.0005);
double a = 0.0;
const int maxGdIt = 9000;
const double delta = 0.009;
const double alpha = 0.0000025;
double constTikh = 2 * alpha;

/*signum function*/
double sign(double n)
{
	if (n > 0.0)
		return 1.0;
	else if (n < 0.0)
		return -1.0;
	else if (n == 0.0)
		return 0.0;
}

class OpticalFlowMethod
{
public:
	OpticalFlowMethod(int sizeX, int sizeY, double t, double dt);
	int GetDegreesOfFreedom();
	void SetW(int& sizeX);
	void SetH(int& sizeY);
	void SetK(int k);
	void SetHX(double _hx);
	void SetHY(double _hy);
	void SetT(double t);
	void SetDT(double _dt);
	int GetW() const;
	int GetH() const;
	int GetK() const;
	double GetHX() const;
	double GetHY() const;
	double GetT() const;
	double GetDT() const;

	void PrintFunc1D(float* f, const char* filename, int gradDes);
	void PrintFunc2D(float* f1, float* f2, const char* filename, int gradDes);

	void SetInitGamma(float* gamma1, float* u, float* I, int gradDes);
	void SetInitFunction(float* f, int gradDes, const char* function, const char* variable);
	void SetInitFunction2(float* f, int gradDes, const char* function, const char* variable);
	void SetInitData(float* f, float* data);
	void SetInitA(float* a1, float* a2);

	void PrimaryU(float* u, float* a1, float* a2, float* data, const char* filename, int gradDes);
	void AdjointGamma(float* u, float* gamma1, float* a1, float* a2, float* I, int gradDes);
	void VelocityUpdate(float* u, float* gamma1, float* a1, float* a2, const double& delta, const double& aplha, int gradDes);

	bool isCFL(float* a1, float* a2, float* cfl1, float* cfl2, int gradDes);

	~OpticalFlowMethod() {}

private:
	int W;			//image width
	int H;			//image height
	double hx;		//step on the x axis
	double hy;		//step on the y axis 

	double T, dt;	//real time, time step
	int K;			//Number of time steps
};

OpticalFlowMethod::OpticalFlowMethod(int sizeX, int sizeY, double t, double dt)
{
	SetW(sizeX);
	SetH(sizeY);
	SetHX(1.0 / (double)W);
	SetHY(1.0 / (double)H);
	SetT(t);
	SetDT(dt);
	SetK(ceil(T / dt));
}

void OpticalFlowMethod::SetW(int& sizeX)
{
	this->W = sizeX;
}
void OpticalFlowMethod::SetH(int& sizeY)
{
	this->H = sizeY;
}
void OpticalFlowMethod::SetK(int k)
{
	this->K = k;
}
void OpticalFlowMethod::SetHX(double _hx)
{
	this->hx = _hx;
}
void OpticalFlowMethod::SetHY(double _hy)
{
	this->hy = _hy;
}
void OpticalFlowMethod::SetT(double t)
{
	this->T = t;
}
void OpticalFlowMethod::SetDT(double _dt)
{
	this->dt = _dt;
}
int OpticalFlowMethod::GetW() const
{
	return this->W;
}
int OpticalFlowMethod::GetH() const
{
	return this->H;
}
int OpticalFlowMethod::GetK() const
{
	return this->K;
}
double OpticalFlowMethod::GetHX() const
{
	return this->hx;
}
double OpticalFlowMethod::GetHY() const
{
	return this->hy;
}
double OpticalFlowMethod::GetT() const
{
	return this->T;
}
double OpticalFlowMethod::GetDT() const
{
	return this->dt;
}

int OpticalFlowMethod::GetDegreesOfFreedom()
{
	return this->W + this->H * this->W + this->W * this->H * (this->K + 1);
}

//Save the scalar data
void OpticalFlowMethod::PrintFunc1D(float* f, const char* filename, int gradDes)
{
	for (int k = 0; k < K; k++)
	{
		if (k % 200 == 0 || k == 499)
		{
			stringstream st;
			st << filename << "-" << k << "_" << gradDes << ".txt";
			fstream file1;

			file1.open(st.str(), ios::out);

			int kWH = k * W * H;
			for (int j = 0; j < H; j++)
			{
				int jW = j * W;
				double jhy = j * hy;
				for (int i = 0; i < W; i++)
				{
					int s = i + jW + kWH;
					file1 << i * hx << "\t" << jhy << "\t \t" << f[s] << endl;
				}
				file1 << endl;
			}
			file1.close();
		}
	}
}

//Save the vector data
void OpticalFlowMethod::PrintFunc2D(float* f1, float* f2, const char* filename, int gradDes)
{
	for (int k = 0; k < K; k++)
	{
		if (k % 200 == 0 || k == 999)
		{
			stringstream st;
			st << filename << "-" << k << "_" << gradDes << ".txt";
			fstream file1;

			file1.open(st.str(), ios::out);

			int kWH = k * W * H;
			for (int j = 0; j < H; j++)
			{
				int jW = j * W;
				double jhy = j * hy;
				for (int i = 0; i < W; i++)
				{
					int s = i + jW + kWH;
					file1 << i * hx << "\t" << jhy << "\t" << f1[s] << "\t \t" << f2[s] << endl;
				}
				file1 << endl;
			}
			file1.close();
		}
	}
}

void OpticalFlowMethod::SetInitFunction2(float* f, int gradDes, const char* function, const char* variable)
{
	double z = 0.35 + (0.0 - 0.5) * 5 * dt / 2.0;
	for (int j = 0; j < H; j++)
	{
		int jW = j * W;
		double y = j * hy;
		for (int i = 0; i < W; i++)
		{
			int s = i + jW;
			double x = i * hx;
			if (i == 0 || j == 0 || i == (W - 1) || j == (H - 1))
			{
				f[s] = 0.0;
				continue;
			}
			f[s] = exp((-((x - z) * (x - z)) - ((y - 0.5) * (y - 0.5))) / 0.03);
		}
	}
}

/*Initial condition - Gauss or signum function*/
void OpticalFlowMethod::SetInitFunction(float* f, int gradDes, const char* function, const char* variable)
{
	for (int k = 0; k < K; k++)
	{
		double z = 0.35 + (k - 0.5) * 5 * dt / 2.0;
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			double y = j * hy;
			for (int i = 0; i < W; i++)
			{
				int s = i + jW + kWH;
				double x = i * hx;
				if (i == 0 || j == 0 || i == (W - 1) || j == (H - 1))
				{
					f[s] = 0.0;
					continue;
				}
				else
				{
					if (function == "gauss")
					{
						f[s] = exp((-((x - z) * (x - z)) - ((y - 0.5) * (y - 0.5))) / 0.03);
					}
					else if (function == "signum")
					{
						double v = (x - 10.0) * (x - 10.0) + (y - 10.0) * (y - 10.0) - 10.0;
						f[s] = sign(v) / 2.0 - 0.5;
					}
					else
					{
						cout << "Function not available for initialization." << endl;
					}
				}
			}
		}
	}
	if (gradDes == 1)
	{
		if (variable == "u")
			PrintFunc1D(f, "outputu0", gradDes);
		else if (variable == "I")
			PrintFunc1D(f, "outputI", gradDes);
	}
}

/*Initial condition*/
void OpticalFlowMethod::SetInitData(float* f, float* data)
{
	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			for (int i = 0; i < W; i++)
			{
				int s = i + jW + kWH;
				f[s] = data[s];
			}
		}
	}
}

/*Initial condition - adjoint problem*/
void OpticalFlowMethod::SetInitGamma(float* gamma1, float* u, float* I, int gradDes)
{
	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			for (int i = 0; i < W; i++)
			{
				int r = i + jW + kWH;
				gamma1[r] = 0.0;
			}
		}
	}
}

/*Initial estimate of optical flow*/
void OpticalFlowMethod::SetInitA(float* a1, float* a2)
{
	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			for (int i = 0; i < W; i++)
			{
				int s = i + jW + kWH;
				a1[s] = a;
				a2[s] = a;
			}
		}
	}
}


/*Primary problem*/
void OpticalFlowMethod::PrimaryU(float* u, float* a1, float* a2, float* data, const char* filename, int gradDes)
{
	SetInitFunction(u, gradDes, "gauss", "u");

	int WH = W * H;
	double dt_hx = dt / hx;
	double dt_hy = dt / hy;

	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			for (int i = 0; i < W; i++)
			{
				int s = i + jW + kWH;
				int sWH = s + WH;

				if (i == 0 || j == 0 || i == (W - 1) || j == (H - 1) || k == 0)
				{
					u[sWH] = u[s];
					continue;
				}
				else
				{
					float us = u[s];
					float a1s = a1[s];
					float a2s = a2[s];

					if (a1s < 0 && a2s < 0)
					{
						u[sWH] = us - (a1s * (dt_hx)) * (u[s + 1] - us) - (a2s * (dt_hy)) * (u[s + W] - us);
					}
					else if (a1s >= 0 && a2s < 0)
					{
						u[sWH] = us - (a1s * (dt_hx)) * (us - u[s - 1]) - (a2s * (dt_hy)) * (u[s + W] - us);
					}
					else if (a1s < 0 && a2s >= 0)
					{
						u[sWH] = us - (a1s * (dt_hx)) * (u[s + 1] - us) - (a2s * (dt_hy)) * (us - u[s - W]);
					}
					else
					{
						u[sWH] = us - (a1s * (dt_hx)) * (us - u[s - 1]) - (a2s * (dt_hy)) * (us - u[s - W]);
					}
				}
			}
		}
	}
	if (gradDes % 1000 == 0 || gradDes == 1)
	{
		PrintFunc1D(u, filename, gradDes);
	}
}

/*Adjoint problem*/
void OpticalFlowMethod::AdjointGamma(float* u, float* gamma1, float* a1, float* a2, float* I, int gradDes)
{
	SetInitGamma(gamma1, u, I, gradDes);

	int WH = W * H;
	double dt_hx = dt / hx;
	double dt_hy = dt / hy;
	double dt2 = dt * 2.0;

	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			for (int i = 0; i < W; i++)
			{
				int r = i + jW + kWH;
				int rWH = r + WH;

				if (i == 0 || j == 0 || i == (W - 1) || j == (H - 1))
				{
					gamma1[rWH] = 0.0;
					continue;
				}
				else
				{
					float gam1r = gamma1[r];
					float a1r = a1[r];
					float a2r = a2[r];

					if (a1r < 0.0 && a2r < 0.0)
					{
						gamma1[rWH] = gam1r - (a1r * (dt_hx)) * (gamma1[r + 1] - gam1r) - (a2r * (dt_hy)) * (gamma1[r + W] - gam1r)
							+ dt2 * (I[r] - u[r]);
					}
					else if (a1r >= 0.0 && a2r < 0.0)
					{
						gamma1[rWH] = gam1r - (a1r * (dt_hx)) * (gam1r - gamma1[r - 1]) - (a2r * (dt_hy)) * (gamma1[r + W] - gam1r)
							+ dt2 * (I[r] - u[r]);
					}
					else if (a1r < 0.0 && a2r >= 0.0)
					{
						gamma1[rWH] = gam1r - (a1r * (dt_hx)) * (gamma1[r + 1] - gam1r) - (a2r * (dt_hy)) * (gam1r - gamma1[r - W])
							+ dt2 * (I[r] - u[r]);
					}
					else
					{
						gamma1[rWH] = gam1r - (a1r * (dt_hx)) * (gam1r - gamma1[r - 1]) - (a2r * (dt_hy)) * (gam1r - gamma1[r - W])
							+ dt2 * (I[r] - u[r]);
					}
				}
			}
		}
	}
	if (gradDes % 1000 == 0)
	{
		PrintFunc1D(gamma1, "outputg", gradDes);
	}
}

/*Updating the optical flow - gradient descent*/
void OpticalFlowMethod::VelocityUpdate(float* u, float* gamma1, float* a1, float* a2, const double& delta, const double& alph, int gradDes)
{
	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			for (int i = 0; i < W; i++)
			{
				int s = i + jW + kWH;
				a1[s] = a1[s] - delta * gamma1[s] * ((u[s + 1] - u[s]) / hx) - constTikh * a1[s];
				a2[s] = a2[s] - delta * gamma1[s] * ((u[s + W] - u[s]) / hy) - constTikh * a2[s];
			}
		}
	}
	if (gradDes % 1000 == 0)
	{
		PrintFunc2D(a1, a2, "outputa", gradDes);
	}
}

bool OpticalFlowMethod::isCFL(float* a1, float* a2, float* cfl1, float* cfl2, int gradDes)
{
	bool b = true;
	double dt_hx = dt / hx;
	double dt_hy = dt / hy;

	for (int k = 0; k < K; k++)
	{
		int kWH = k * W * H;
		for (int j = 0; j < H; j++)
		{
			int jW = j * W;
			double jhy = j * hy;
			for (int i = 0; i < W; i++)
			{
				int s = i + jW + kWH;
				cfl1[s] = a1[s] * (dt_hx);
				cfl2[s] = a2[s] * (dt_hy);
				if (abs(cfl1[s]) > 1.0 || abs(cfl2[s]) > 1.0)
				{
					cout << "CFL porusena v k = " << k << " na i * hx: " << i * hx << ", j * hy: " << jhy << ", cfl1[s] je: " << cfl1[s] << ", cfl2[s] je: " << cfl2[s] << endl;
					b = false;
				}
			}
		}
	}
	if (gradDes % 1000 == 0)
	{
		PrintFunc2D(cfl1, cfl2, "outputcfl", gradDes);
	}
	if (b == false)
		return false;
	else
		return true;
}

int main()
{
	system("CHCP 1250 > NUL");
	auto start = high_resolution_clock::now();
	OpticalFlowMethod problem(sizes, sizes, finalTime, timeStep);

	int deg = problem.GetDegreesOfFreedom();
	//cout << deg << endl;

	float* u = new float[deg];
	float* gamma1 = new float[deg];
	float* a1 = new float[deg];
	float* a2 = new float[deg];
	float* I = new float[deg];
	float* cfl1 = new float[deg];
	float* cfl2 = new float[deg];

	/*Initialize the velocity vector*/
	problem.SetInitA(a1, a2);
	problem.SetInitFunction(I, 1, "gauss", "I");

	/*Gradient descent iteration*/
	for (int gdIt = 1; gdIt <= maxGdIt; gdIt++)
	{
		problem.PrimaryU(u, a1, a2, I, "outputu", gdIt);
		problem.AdjointGamma(u, gamma1, a1, a2, I, gdIt);
		problem.VelocityUpdate(u, gamma1, a1, a2, delta, alpha, gdIt);
		if (!problem.isCFL(a1, a2, cfl1, cfl2, gdIt))
			break;

		if (gdIt % 10 == 0)
		{
			cout << "Konec iterace èíslo " << gdIt << endl;
			cout << "*______________________*" << endl;
		}
		if (gdIt % 100 == 0 || gdIt == 1)
		{
			auto stopp = high_resolution_clock::now();
			auto during = duration_cast<microseconds>(stopp - start);
			auto tMinutes = during / 60000000;
			auto tSeconds = during / 1000000;
			cout << "The runtime in the " << gdIt << ". iteration is " << tMinutes.count() << " minutes or " << tSeconds.count() << " seconds." << endl;
		}
	}

	auto stop = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(stop - start);
	auto timeMinutes = duration / 60000000;
	auto timeSeconds = duration / 1000000;
	cout << "The runtime is " << timeMinutes.count() << " minutes or " << timeSeconds.count() << " seconds." << endl;
}