#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
using namespace std;


inline const char* const BoolToString(bool& b)
{
	return b ? "true" : "false";
}
class Node {
public:

	double x;
	double y;
	bool BC;
	Node() :x(0), y(0)
	{
	};
	Node(double _x, double _y) : x(_x), y(_y)
	{

	};
	Node(double _x, double _y, bool _BC) : x(_x), y(_y), BC(_BC)
	{


	};
	Node& operator=(const Node& other)
	{
		this->x = other.x;
		this->y = other.y;
		this->BC = other.BC;
		return *this;
	}

};

class Element
{
public:
	int id[4];
	double k = 25.0;
	double alfa = 300.0;
	double ambient_temp = 1200.0;
	double c = 700.0;
	double ro = 7800;
	double dtau = 50;
	//bool bc;
	vector<vector<double>> matrixH;
	vector<vector<double>> matrixHBC;
	vector<vector<double>> matrixC;
	vector<double> vectorP;
	void printMatrixH()
	{
		for (int i = 0; i < this->matrixH.size(); i++)
		{
			for (int j = 0; j < this->matrixH[i].size(); j++)
			{
				cout << this->matrixH[i][j] << " ";

			}
			cout << endl;
		}
	}
	void printMatrixHBC()
	{
		for (int i = 0; i < this->matrixHBC.size(); i++)
		{
			for (int j = 0; j < this->matrixHBC[i].size(); j++)
			{
				cout << this->matrixHBC[i][j] << " ";
			}
			cout << endl;
		}
	}
	void printMatrixC()
	{
		for (int i = 0; i < this->matrixC.size(); i++)
		{
			for (int j = 0; j < this->matrixC[i].size(); j++)
			{
				cout << this->matrixC[i][j] << " ";
			}
			cout << endl;
		}
	}
	void printVectorP()
	{
		for (int i = 0; i < this->vectorP.size(); i++)
		{

			cout << this->vectorP[i] << " ";


		}
		cout << endl;
	}

};

class Grid
{
	double H;
	double B;

	int nN;
	int nE;


public:
	int nH;
	int nB;
	vector<vector<Node>> nodes;
	vector <Element> elements;
	Grid(double h, double b, int nh, int nb) : H(h), B(b), nH(nh), nB(nb)
	{
		nN = nH * nB;
		nE = (nH - 1) * (nB - 1);
		double dh = H / (nH - 1);
		double db = B / (nB - 1);
		elements.resize(nE);
		nodes.resize(nH);

		for (int i = 0; i < nH; i++)
			nodes[i].resize(nB);

		cout << "dh = " << dh << " db: " << db << endl;

		for (int i = 0; i < nH; i++)
		{
			for (int j = 0; j < nB; j++)
			{

				if (i == 0 || i == nH - 1 || j == 0 || j == nB - 1)
					nodes[i][j] = Node(j * db, i * dh, true);
				else
					nodes[i][j] = Node(j * db, i * dh, false);

			}
		}

		int helper = 0;
		for (int i = 0; i < nE; i++)
		{
			if (i != 0)
				if (i % (nH - 1) == 0)
					helper++;
			elements[i].id[0] = i + helper;
			elements[i].id[1] = elements[i].id[0] + nH;
			elements[i].id[2] = elements[i].id[1] + 1;
			elements[i].id[3] = elements[i].id[0] + 1;


		}

	}

	Node getNodeFromElement(int elementID, int tabID)
	{

		return getNode(elements[elementID].id[tabID]);
	}
	Node getNode(int nodeID)
	{
		int j = nodeID / nH;
		int i = nodeID % nH;


		return nodes[i][j];
	}
	int getNumberOfelements()
	{
		return elements.size();
	}
	int getNumberOfNodes()
	{
		return nH * nB;
	}
	void print()
	{

		cout << "Printing elements and their nodes:\n";
		for (int i = 0; i < nE; i++)
		{
			cout << "Element " << i + 1 << " nodes: " << elements[i].id[0] + 1 << ", " << elements[i].id[1] + 1 << ", " << elements[i].id[2] + 1 << ", " << elements[i].id[3] + 1 << endl;
		}
		cout << "Printing nodes and their position:\n";

		for (int j = 0; j < nB;j++)
		{
			for (int i = 0; i < nH; i++)
			{
				int id = nH * j + i + 1;
				cout << "Node " << id << " x: " << nodes[i][j].x << ", y: " << nodes[i][j].y << ", BC: " << BoolToString(nodes[i][j].BC) << endl;
			}
		}
	}
};
class Wagi {
public:

	vector<double> wezly;
	vector<double> wspolczynniki;

};

class Wagi2pukty : public Wagi {
public:

	Wagi2pukty()
	{
		wezly.push_back(-1 / sqrt(3));
		wezly.push_back(1 / sqrt(3));
		wspolczynniki.push_back(1);
		wspolczynniki.push_back(1);
	}

};
class Wagi3pukty : public Wagi {
public:


	Wagi3pukty()
	{

		wezly.push_back(-sqrt(3.0 / 5.0));
		wezly.push_back(0);
		wezly.push_back(sqrt(3.0 / 5.0));
		wspolczynniki.push_back(5.0 / 9.0);
		wspolczynniki.push_back(8.0 / 9.0);
		wspolczynniki.push_back(5.0 / 9.0);
	}

};
class Element4
{
public:
	vector<vector<double>> dKsi;
	vector<vector<double>> dNi;


	Element4(Wagi punkty)
	{

		dKsi.resize(punkty.wezly.size() * punkty.wezly.size());
		dNi.resize(punkty.wezly.size() * punkty.wezly.size());
		for (int i = 0; i < dKsi.size(); i++)
		{
			dKsi[i].resize(4);
			dNi[i].resize(4);
		}

		for (int i = 0; i < punkty.wezly.size(); i++)
		{

			for (int j = 0; j < punkty.wezly.size(); j++)
			{
				dKsi[i * punkty.wezly.size() + j][0] = -1 / 4.0 * (1 - punkty.wezly[i]) * punkty.wspolczynniki[i];
				dKsi[i * punkty.wezly.size() + j][1] = 1 / 4.0 * (1 - punkty.wezly[i]) * punkty.wspolczynniki[i];
				dKsi[i * punkty.wezly.size() + j][2] = 1 / 4.0 * (1 + punkty.wezly[i]) * punkty.wspolczynniki[i];
				dKsi[i * punkty.wezly.size() + j][3] = -1 / 4.0 * (1 + punkty.wezly[i]) * punkty.wspolczynniki[i];
				dNi[i * punkty.wezly.size() + j][0] = -1 / 4.0 * (1 - punkty.wezly[j]) * punkty.wspolczynniki[j];
				dNi[i * punkty.wezly.size() + j][1] = -1 / 4.0 * (1 + punkty.wezly[j]) * punkty.wspolczynniki[j];
				dNi[i * punkty.wezly.size() + j][2] = 1 / 4.0 * (1 + punkty.wezly[j]) * punkty.wspolczynniki[j];
				dNi[i * punkty.wezly.size() + j][3] = 1 / 4.0 * (1 - punkty.wezly[j]) * punkty.wspolczynniki[j];
			}

		}
	}

	void print()
	{
		cout << "Ksi: " << endl;
		for (int i = 0; i < dKsi.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << dKsi[i][j] << " ";
			}
			cout << endl;

		}

		cout << "Ni: " << endl;
		for (int i = 0; i < dNi.size(); i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << dNi[i][j] << " ";
			}
			cout << endl;

		}
	}
};
vector<vector<double>> jacobian(int idElementu, int j, Grid& grid, Element4 element, double& detJ)
{
	int jacobianSize = 2;
	vector<vector<double>> jakobianMatrix;
	jakobianMatrix.resize(jacobianSize);
	vector<vector<double>> inverseJakobianMatrix;
	inverseJakobianMatrix.resize(jacobianSize);
	for (int k = 0; k < jacobianSize; k++)
	{
		jakobianMatrix[k].resize(jacobianSize);
		inverseJakobianMatrix[k].resize(jacobianSize);
	}

	for (int k = 0; k < 4; k++)
	{
		jakobianMatrix[0][0] += element.dKsi[j][k] * grid.getNodeFromElement(idElementu, k).x;
		jakobianMatrix[0][1] += element.dKsi[j][k] * grid.getNodeFromElement(idElementu, k).y;
		jakobianMatrix[1][0] += element.dNi[j][k] * grid.getNodeFromElement(idElementu, k).x;
		jakobianMatrix[1][1] += element.dNi[j][k] * grid.getNodeFromElement(idElementu, k).y;
	}

	double detJacobian = (jakobianMatrix[0][0] * jakobianMatrix[1][1]) - (jakobianMatrix[0][1] * jakobianMatrix[1][0]);
	detJ = detJacobian;
	inverseJakobianMatrix[0][0] = jakobianMatrix[1][1] / detJacobian;
	inverseJakobianMatrix[0][1] = jakobianMatrix[0][1] / detJacobian;
	inverseJakobianMatrix[0][1] *= -1;
	/*if (0 == inverseJakobianMatrix[0][1])
		inverseJakobianMatrix[0][1] = 0;
	inverseJakobianMatrix[1][0] = jakobianMatrix[1][0] / detJacobian;
	inverseJakobianMatrix[1][0] *= -1;
	if (0 == inverseJakobianMatrix[1][0])
		inverseJakobianMatrix[1][0] = 0;*/
	inverseJakobianMatrix[1][1] = jakobianMatrix[0][0] / detJacobian;
	/*cout << "jacobian matrix : " << endl;
	for (int i = 0; i < jakobianMatrix.size(); i++)
	{
		for (int j = 0; j < jakobianMatrix.size(); j++)
			cout << jakobianMatrix[i][j] << " ";
		cout << endl;
	}*/

	return inverseJakobianMatrix;

}

void Gauss(Wagi  a, int wymiar)
{
	double wynik = 0;
	double wynik2 = 0;
	if (wymiar == 1)
	{
		for (int i = 0; i < a.wspolczynniki.size(); i++)
		{
			wynik += a.wspolczynniki[i] * (5 * (a.wezly[i]) * (a.wezly[i]) + 3 * a.wezly[i] + 6);
		}
		cout.precision(15);
		cout << "Wynik calki dla przestrzeni " << wymiar << "d stosujac " << a.wezly.size() << " punktowy schemat calkowania to: " << wynik << endl;
	}
	else if (wymiar == 2)
	{
		for (int i = 0; i < a.wspolczynniki.size(); i++)
		{
			for (int j = 0; j < a.wspolczynniki.size(); j++)
			{
				wynik2 += a.wspolczynniki[i] * a.wspolczynniki[j] * (5 * (a.wezly[i]) * (a.wezly[i]) * (a.wezly[j]) * (a.wezly[j]) + 3 * a.wezly[i] * a.wezly[j] + 6);
			}
		}
		cout.precision(15);
		cout << "Wynik calki dla przestrzeni " << wymiar << "d stosujac " << a.wezly.size() << " punktowy schemat calkowania to: " << wynik2 << endl;
	}
	else {
		cout << "Z�y schemat \n";
	}
	return;

}
void dNdxXdNdY(vector<vector<double>>& dNdX, vector<vector<double>>& dNdY, vector<vector<double>>& inversejacobian, Element4& c, int j)
{
	dNdX[j][0] = inversejacobian[0][0] * c.dKsi[j][0] + inversejacobian[0][1] * c.dNi[j][0];
	dNdX[j][1] = inversejacobian[0][0] * c.dKsi[j][1] + inversejacobian[0][1] * c.dNi[j][1];
	dNdX[j][2] = inversejacobian[0][0] * c.dKsi[j][2] + inversejacobian[0][1] * c.dNi[j][2];
	dNdX[j][3] = inversejacobian[0][0] * c.dKsi[j][3] + inversejacobian[0][1] * c.dNi[j][3];
	dNdY[j][0] = inversejacobian[1][0] * c.dKsi[j][0] + inversejacobian[1][1] * c.dNi[j][0];
	dNdY[j][1] = inversejacobian[1][0] * c.dKsi[j][1] + inversejacobian[1][1] * c.dNi[j][1];
	dNdY[j][2] = inversejacobian[1][0] * c.dKsi[j][2] + inversejacobian[1][1] * c.dNi[j][2];
	dNdY[j][3] = inversejacobian[1][0] * c.dKsi[j][3] + inversejacobian[1][1] * c.dNi[j][3];
}

void calculatePieceOfHmatrix(vector<double>& dNdX, vector<double>& dNdY, vector<vector<double>>& H, double detJ, double k)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			double tempX = dNdX[i] * dNdX[j];
			double tempY = dNdY[i] * dNdY[j];
			H[i][j] = (tempX + tempY) * detJ * k;

		}
	}

}vector<pair<double, double>> makePoints(Element4 el, vector<pair<double, double>>& wsp) {
	Wagi w;
	if (el.dKsi.size() > 4)
		w = Wagi3pukty();
	else
		w = Wagi2pukty();

	vector<pair<double, double>> edgePoints;
	edgePoints.resize(w.wezly.size() * w.wezly.size());
	wsp.resize(w.wezly.size() * w.wezly.size());

	for (int i = 0; i < w.wezly.size(); i++) {
		for (int j = 0; j < w.wezly.size(); j++) {
			edgePoints[i * w.wezly.size() + j].first = w.wezly[i];
			edgePoints[i * w.wezly.size() + j].second = w.wezly[j];
			wsp[i * w.wezly.size() + j].first = w.wspolczynniki[i];
			wsp[i * w.wezly.size() + j].second = w.wspolczynniki[j];
			cout << "ksi: " << edgePoints[i * w.wezly.size() + j].first << " ,Ni: " << edgePoints[i * w.wezly.size() + j].second << " , waga1: " << wsp[i * w.wezly.size() + j].first << " , waga2: " << wsp[i * w.wezly.size() + j].second << endl;
		}
	}
	return edgePoints;
}
vector<double>	calculatShapeFunctions(pair<double, double>& edgePoints) {

	vector<double> shapeFunctions;
	shapeFunctions.resize(4);
	shapeFunctions[0] = 0.25 * (1.0 - edgePoints.first) * (1.0 - edgePoints.second);
	shapeFunctions[1] = 0.25 * (1.0 + edgePoints.first) * (1.0 - edgePoints.second);
	shapeFunctions[2] = 0.25 * (1.0 + edgePoints.first) * (1.0 + edgePoints.second);
	shapeFunctions[3] = 0.25 * (1.0 - edgePoints.first) * (1.0 + edgePoints.second);

	return shapeFunctions;
}

void calculatePeaceOfPieceOfC(vector<vector<double>>& H, vector<double>& shapeFunctions, pair<double, double>& wages, double detJ, double c, double ro) {
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{


			H[i][j] = (shapeFunctions[i] * shapeFunctions[j]) * detJ * c * ro/** wages.first*//**wages.second*/;

		}
	}
}
void makeSquareMatrixOfSize(vector<vector<double>>& H, int n)
{
	H.resize(n);
	for (int i = 0; i < H.size(); i++) {
		H[i].resize(n);
		for (int j = 0; j < H[i].size(); j++)
		{
			H[i][j] = 0.0;
		}
	}
}
void addMatrixH(vector<vector<double>>& H1, vector<vector<double>>& H2)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			H2[i][j] += H1[i][j];
	}
}
void zeroMatrixH(vector<vector<double>>& H)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			H[i][j] = 0;
	}
}
void resizeMatrixesForElementSize(vector<vector<double>>& H1, vector<vector<double>>& H2, vector<vector<double>>& dNdX, vector<vector<double>>& dNdY, vector<vector<double>>& Hbc1, vector<vector<double>>& Hbc2, Element4& c)
{
	H1.resize(4);
	H2.resize(4);
	Hbc1.resize(4);
	Hbc2.resize(4);
	dNdX.resize(c.dNi.size());
	dNdY.resize(c.dNi.size());
	for (int k = 0; k < c.dNi.size(); k++)
	{
		dNdY[k].resize(4);
		dNdX[k].resize(4);

	}
	for (int k = 0; k < 4; k++)
	{
		H1[k].resize(4);
		H2[k].resize(4);
		Hbc1[k].resize(4);
		Hbc2[k].resize(4);
		for (int j = 0; j < 4; j++)
		{
			H2[k][j] = 0.0;
			Hbc2[k][j] = 0.0;
		}
	}
}
vector<vector<pair<double, double>>> makeEdgePoints(Wagi w) {


	vector<vector<pair<double, double>>> edgePoints;
	edgePoints.resize(4);
	for (int i = 0; i < 4; i++)
	{
		edgePoints[i].resize(w.wezly.size());
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < w.wezly.size(); j++) {
			if (i == 0)
			{
				edgePoints[i][j].first = -1.0;
				edgePoints[i][j].second = w.wezly[j];

			}
			else if (i == 1)
			{
				edgePoints[i][j].second = -1.0;
				edgePoints[i][j].first = w.wezly[j];

			}
			else if (i == 2)
			{
				edgePoints[i][j].first = 1.0;
				edgePoints[i][j].second = w.wezly[j];

			}
			else
			{
				edgePoints[i][j].second = 1.0;
				edgePoints[i][j].first = w.wezly[j];

			}
		}
	}
	return edgePoints;
}

vector<double>	calculateEdgeShapeFunctions(pair<double, double>& edgePoints) {

	vector<double> shapeFunctions;
	shapeFunctions.resize(4);
	shapeFunctions[0] = 0.25 * (1.0 - edgePoints.first) * (1.0 - edgePoints.second);
	shapeFunctions[1] = 0.25 * (1.0 + edgePoints.first) * (1.0 - edgePoints.second);
	shapeFunctions[2] = 0.25 * (1.0 + edgePoints.first) * (1.0 + edgePoints.second);
	shapeFunctions[3] = 0.25 * (1.0 - edgePoints.first) * (1.0 + edgePoints.second);

	return shapeFunctions;
}
void calculatePeaceOfPieceOfHbc(vector<vector<double>>& H, vector<double>& shapeFunctions, double weight, double detJ, double alfa) {
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{


			H[i][j] = (shapeFunctions[i] * shapeFunctions[j]) * weight * detJ * alfa;

		}
	}

}
void zeroVector(vector<double>& p) {
	for (int i = 0; i < p.size(); i++) {
		p[i] = 0.0;
	}
}
void makeVectorOfsize(vector<double>& p, int n) {
	p.resize(n);
	for (int i = 0; i < p.size(); i++) {
		p[i] = 0.0;
	}
}

double calculateDetJForHBC(Grid& grid, int index, int i)
{
	double detJ = 0.1;
	switch (i)
	{
	case 0:
		if (grid.getNodeFromElement(index, 0).BC == true && grid.getNodeFromElement(index, 3).BC == true)
			detJ = abs(grid.getNodeFromElement(index, 0).y - grid.getNodeFromElement(index, 3).y);
		else
			detJ = 0;
		break;
	case 1:
		if (grid.getNodeFromElement(index, 0).BC == true && grid.getNodeFromElement(index, 1).BC == true)
			detJ = abs(grid.getNodeFromElement(index, 0).x - grid.getNodeFromElement(index, 1).x);
		else
			detJ = 0;
		break;
	case 2:
		if (grid.getNodeFromElement(index, 1).BC == true && grid.getNodeFromElement(index, 2).BC == true)
			detJ = abs(grid.getNodeFromElement(index, 1).y - grid.getNodeFromElement(index, 2).y);
		else
			detJ = 0;
		break;
	case 3:
		if (grid.getNodeFromElement(index, 2).BC == true && grid.getNodeFromElement(index, 3).BC == true)
			detJ = abs(grid.getNodeFromElement(index, 2).x - grid.getNodeFromElement(index, 3).x);
		else
			detJ = 0;
		break;

	default:
		break;
	}
	return detJ;

}

void calculateHbc(vector<vector<double>>& Hbc1, vector<vector<double>>& Hbc2, Element4 el, double alfa, Grid& grid, int index) {

	vector<vector<double>> Hbc3;
	makeSquareMatrixOfSize(Hbc3, 4);

	Wagi c;
	if (el.dKsi.size() > 4)
		c = Wagi3pukty();
	else
		c = Wagi2pukty();
	zeroMatrixH(Hbc1);
	zeroMatrixH(Hbc2);
	zeroMatrixH(Hbc3);
	vector<vector<pair<double, double>>> edgePoints = makeEdgePoints(c);
	for (int i = 0; i < edgePoints.size(); i++) {
		for (int j = 0; j < edgePoints[i].size(); j++) {


			vector<double> shapeFunctions = calculateEdgeShapeFunctions(edgePoints[i][j]);
			double detJ = calculateDetJForHBC(grid, index, i);
			detJ /= 2.0;
			//cout << "ksi: " << edgePoints[i][j].first << " ni: " << edgePoints[i][j].second << " N1: "<<shapeFunctions[0] << " N2: " << shapeFunctions[1]<< " N3: " << shapeFunctions[2]<< " N4: " << shapeFunctions[3]<< " Alfa: "<<alfa<< " waga: " << c.wspolczynniki[j]<<" detJ "<< detJ <<endl;

			//cout << "detJ: " << detJ << " ";
			calculatePeaceOfPieceOfHbc(Hbc3, shapeFunctions, c.wspolczynniki[j], detJ, alfa);

			addMatrixH(Hbc3, Hbc1);
			zeroMatrixH(Hbc3);
		}
		addMatrixH(Hbc1, Hbc2);

		zeroMatrixH(Hbc1);

	}

}


void calculatePeaceOfPieceOfVectorP(vector<double>& H, vector<double>& shapeFunctions, double weight, double detJ, double alfa, double ambient_temp) {
	for (int i = 0; i < 4; i++)
	{


		H[i] = (shapeFunctions[i]) * weight * detJ * alfa * ambient_temp;


	}

}
void addVectors(vector<double>& p1, vector<double>& p2) {
	for (int i = 0; i < p1.size(); i++) {
		p2[i] += p1[i];
	}
}


void calculateVectorP(vector<double>& p1, vector<double>& p2, Element4 el, double alfa, Grid& grid, int index, double ambient_temp) {

	vector<double> p3;
	makeVectorOfsize(p3, 4);

	Wagi c;
	if (el.dKsi.size() > 4)
		c = Wagi3pukty();
	else
		c = Wagi2pukty();
	zeroVector(p1);
	zeroVector(p2);
	zeroVector(p3);
	vector<vector<pair<double, double>>> edgePoints = makeEdgePoints(c);
	for (int i = 0; i < edgePoints.size(); i++) {
		for (int j = 0; j < edgePoints[i].size(); j++) {


			vector<double> shapeFunctions = calculateEdgeShapeFunctions(edgePoints[i][j]);
			double detJ = calculateDetJForHBC(grid, index, i);
			detJ /= 2.0;
			//	cout << "ksi: " << edgePoints[i][j].first << " ni: " << edgePoints[i][j].second << " N1: " << shapeFunctions[0] << " N2: " << shapeFunctions[1] << " N3: " << shapeFunctions[2] << " N4: " << shapeFunctions[3] << " Alfa: " << alfa << " waga: " << c.wspolczynniki[j] << " detJ " << detJ << endl;

				//cout << "detJ: " << detJ << " ";
			calculatePeaceOfPieceOfVectorP(p3, shapeFunctions, c.wspolczynniki[j], detJ, alfa, ambient_temp);

			addVectors(p3, p1);
			zeroVector(p3);
		}
		addVectors(p1, p2);

		zeroVector(p1);

	}
}
void printVector(vector<double>& s) {
	cout << "vector values: " << endl;
	for (auto x : s)
		cout << x << "  ";
	cout << endl;
}

class ThingsToCount
{

public:
	double dtau = 50;
	double dt = 50;
	double t0 = 100;
	vector<vector<double>> aggregatedHmatrix;
	vector<vector<double>> aggregatedCmatrix;
	vector<vector<double>> aggregatedCmatrixPlusH;
	vector<double> aggregatedVectorP;
	vector<vector<double>> solvedAggregatedHmatrix;
	vector<double> equationSolution;
	vector<double> aggregatedOldVectorP;

	ThingsToCount(int numOfnodes) {
		aggregatedHmatrix.resize(numOfnodes);
		aggregatedCmatrix.resize(numOfnodes);
		solvedAggregatedHmatrix.resize(numOfnodes);
		aggregatedCmatrixPlusH.resize(numOfnodes);
		equationSolution.resize(numOfnodes);                     // here

		for (int i = 0; i < aggregatedHmatrix.size(); i++)
		{
			aggregatedHmatrix[i].resize(numOfnodes);
			aggregatedCmatrix[i].resize(numOfnodes);
			aggregatedCmatrixPlusH[i].resize(numOfnodes);
			solvedAggregatedHmatrix[i].resize(numOfnodes);
			equationSolution[i] = t0;
			for (int j = 0; j < aggregatedHmatrix[i].size(); j++)
			{
				aggregatedHmatrix[i][j] = 0.0;
				aggregatedCmatrix[i][j] = 0.0;
				solvedAggregatedHmatrix[i][j] = 0.0;
				aggregatedCmatrixPlusH[i][j] = 0.0;
			}

		}
		makeVectorOfsize(aggregatedVectorP, numOfnodes);
		makeVectorOfsize(aggregatedOldVectorP, numOfnodes);

	}

	void addHandC() {
		for (int i = 0; i < aggregatedCmatrix.size(); i++) {

			for (int j = 0; j < aggregatedCmatrix[i].size(); j++) {

				aggregatedCmatrixPlusH[i][j] = aggregatedCmatrix[i][j] / dtau + aggregatedHmatrix[i][j];

			}

		}

	}

	void printAggregatedHmatrix() {
		cout << "Aggregated H matrix: " << endl;

		for (int i = 0; i < aggregatedHmatrix.size(); i++) {

			for (int j = 0; j < aggregatedHmatrix[i].size(); j++) {

				cout << aggregatedHmatrix[i][j] << "  ";
			}
			cout << endl;
		}
	}
	void printAggregatedCmatrix() {
		cout << "Aggregated C matrix: " << endl;

		for (int i = 0; i < aggregatedCmatrix.size(); i++) {

			for (int j = 0; j < aggregatedCmatrix[i].size(); j++) {

				cout << aggregatedCmatrix[i][j] << "  ";
			}
			cout << endl;
		}
	}
	void printAggregatedCmatrixPlusH() {
		cout << "Aggregated C + H matrix: " << endl;

		for (int i = 0; i < aggregatedCmatrixPlusH.size(); i++) {

			for (int j = 0; j < aggregatedCmatrixPlusH[i].size(); j++) {

				cout << aggregatedCmatrixPlusH[i][j] << "  ";
			}
			cout << endl;
		}
	}
	void printAggregatedVectorP() {
		cout << "Aggregated vector P: " << endl;

		for (int i = 0; i < aggregatedVectorP.size(); i++) {

			cout << aggregatedVectorP[i] << "  ";
		}
		cout << endl;
	}

	void printWholeSolvedAggregatedHmatrix() {
		cout << "Whole solved aggregated H matrix: " << endl;

		for (int i = 0; i < solvedAggregatedHmatrix.size(); i++) {

			for (int j = 0; j < solvedAggregatedHmatrix[i].size(); j++) {

				cout << solvedAggregatedHmatrix[i][j] << "  ";
			}
			cout << endl;
		}
	}
	void printEquationSolution() {
		cout << "Equation solution: " << endl;

		for (int i = 0; i < equationSolution.size(); i++) {

			cout << "x[" << i << "] = " << equationSolution[i] << "  ";
		}
		cout << endl;
	}
	void vectorPequation()                          //////here
	{
		vector<double> helper;
		helper.resize(aggregatedVectorP.size());
		double temp = 0;
		for (int i = 0; i < aggregatedCmatrix.size(); i++)
		{
			for (int j = 0; j < aggregatedCmatrix[i].size(); j++)
			{

				temp += aggregatedCmatrix[i][j] / dtau * equationSolution[j];

			}
			helper[i] = temp;
			temp = 0;
			aggregatedVectorP[i] = aggregatedOldVectorP[i] + helper[i];
		}
		printVector(helper);
		/*for (int i = 0; i < aggregatedVectorP.size(); i++)
		{
			aggregatedVectorP[i] = aggregatedOldVectorP[i] + helper[i];
		}*/
		//dtau += dt;
	}
	//void solveAggregatedHmatrix() {
	//	solvedAggregatedHmatrix = aggregatedHmatrix;
	//	int n = solvedAggregatedHmatrix.size()+1;
	//	for (int i = 0; i < n - 1; i++)
	//		solvedAggregatedHmatrix[i].push_back(aggregatedVectorP[i]);
	//	for (int i = 0; i < n - 1; i++)
	//	{
	//		for (int j = i + 1; j < n - 1; j++)
	//		{
	//			if (solvedAggregatedHmatrix[i][i] == 0)
	//			{
	//				cout << "Cannot solve the matrix!" << endl;
	//				return ;

	//			}
	//			double helper = (solvedAggregatedHmatrix[j][i] / solvedAggregatedHmatrix[i][i]);
	//			for (int k = 0; k < n; k++)
	//			{
	//				solvedAggregatedHmatrix[j][k] -= helper * solvedAggregatedHmatrix[i][k];
	//			}
	//		}
	//	}
	//	for (int i = n - 2; i >= 0; i--)
	//	{
	//		for (int j = 0; j < n - 1; j++) {
	//			if(i!=j)
	//				solvedAggregatedHmatrix[i][n-1] -= solvedAggregatedHmatrix[i][j] * solvedAggregatedHmatrix[j][n-1];
	//			
	//		}
	//		if (solvedAggregatedHmatrix[i][i] != 0 && solvedAggregatedHmatrix[i][n - 1] != 0)
	//			solvedAggregatedHmatrix[i][n - 1] /= solvedAggregatedHmatrix[i][i];
	//		else if (solvedAggregatedHmatrix[i][i] == 0 && solvedAggregatedHmatrix[i][n - 1] != 0)
	//		{
	//			cout << "Inconsistent system of equations!" << endl;
	//			return;

	//		}
	//	}

	//}


};
//here
void solveAggregatedHmatrix(vector<vector<double>>& aggregatedHmatrix, vector<vector<double>>& solvedAggregatedHmatrix, vector<double>& aggregatedVectorP, vector<double>& equationSolution) {
	int  i, j, k;
	solvedAggregatedHmatrix = aggregatedHmatrix;
	int n = solvedAggregatedHmatrix.size();

	equationSolution.resize(n);
	cout << "OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK LETSSSSSSSSG O                  OOO\n";
	for (i = 0; i < aggregatedVectorP.size(); i++)
		cout << aggregatedVectorP[i] <<" ";


	for (int i = 0; i < n; i++)
		solvedAggregatedHmatrix[i].push_back(aggregatedVectorP[i]);

	for (i = 0; i < n; i++)
		for (k = i + 1; k < n; k++)
			if (abs(solvedAggregatedHmatrix[i][i]) < abs(solvedAggregatedHmatrix[k][i]))
				for (j = 0; j <= n; j++)
				{
					double temp = solvedAggregatedHmatrix[i][j];
					solvedAggregatedHmatrix[i][j] = solvedAggregatedHmatrix[k][j];
					solvedAggregatedHmatrix[k][j] = temp;
				}


	for (i = 0; i < n - 1; i++)
		for (k = i + 1; k < n; k++)
		{
			double t = solvedAggregatedHmatrix[k][i] / solvedAggregatedHmatrix[i][i];
			for (j = 0; j <= n; j++)
				solvedAggregatedHmatrix[k][j] = solvedAggregatedHmatrix[k][j] - t * solvedAggregatedHmatrix[i][j];
		}


	for (i = n - 1; i >= 0; i--)
	{
		equationSolution[i] = solvedAggregatedHmatrix[i][n];
		for (j = i + 1; j < n; j++)
			if (j != i)
				equationSolution[i] = equationSolution[i] - solvedAggregatedHmatrix[i][j] * equationSolution[j];
		equationSolution[i] = equationSolution[i] / solvedAggregatedHmatrix[i][i];
	}

}

void aggregateMatrixes(ThingsToCount& thingsToCount, Grid& grid)
{


	for (auto element : grid.elements) {
		// element.printMatrixH();

		for (int i = 0; i < element.matrixH.size(); i++) {

			for (int j = 0; j < element.matrixH.size(); j++)
			{

				thingsToCount.aggregatedHmatrix[element.id[i]][element.id[j]] += element.matrixH[i][j];
				thingsToCount.aggregatedHmatrix[element.id[i]][element.id[j]] += element.matrixHBC[i][j];
				thingsToCount.aggregatedCmatrix[element.id[i]][element.id[j]] += element.matrixC[i][j];




			}
		}
		//	thingsToCount.printAggregatedHmatrix();
	}


}
void aggregateVectorP(ThingsToCount& thingsToCount, Grid& grid)
{


	for (auto element : grid.elements) {
		// element.printMatrixH();

		for (int i = 0; i < element.vectorP.size(); i++) {


			thingsToCount.aggregatedVectorP[element.id[i]] += element.vectorP[i];
			thingsToCount.aggregatedOldVectorP[element.id[i]] += element.vectorP[i];

		}
		//	thingsToCount.printAggregatedHmatrix();
	}
}

int main()
{
	Grid grid(0.1f, 0.1f, 4, 4);

	//grid.print();
	///*for(int i= 0; i<grid.getNumberOfNodes();i++)
	//	cout << "ID " << i + 1 << " x: " << grid.getNode(i).x << " y: " << grid.getNode(i).y << endl;*/
	//cout << grid.getNodeFromElement(4,0).x <<" , y: " << grid.getNodeFromElement(4,0).y <<endl;
	//cout << grid.getNodeFromElement(4, 1).x << " , y: " << grid.getNodeFromElement(4, 1).y << endl;
	//cout << grid.getNodeFromElement(4, 2).x << " , y: " << grid.getNodeFromElement(4, 2).y << endl;
	//cout << grid.getNodeFromElement(4, 3).x << " , y: " << grid.getNodeFromElement(4, 3).y << endl;
	Wagi  trzyPunkty = Wagi3pukty();
	Wagi  dwaPunkty = Wagi2pukty();

	//Gauss(dwaPunkty, 1);

	Element4 c = Element4(trzyPunkty);




	cout << "dKsi:" << endl;
	for (int i = 0; i < c.dKsi.size(); i++)
	{
		for (int j = 0; j < c.dKsi[i].size(); j++)
		{
			cout << c.dKsi[i][j] << " ";
		}
		cout << endl;

	}
	cout << "dNi:" << endl;
	for (int i = 0; i < c.dNi.size(); i++)
	{
		for (int j = 0; j < c.dNi[i].size(); j++)
		{
			cout << c.dNi[i][j] << " ";
		}
		cout << endl;

	}

	vector<vector<double>> inversejacobian;
	double detJ = 0.0;
	vector<vector<double>> H1;
	vector<vector<double>> H2;
	vector<vector<double>> Hbc1;
	vector<vector<double>> Hbc2;
	vector<vector<double>> dNdX;
	vector<vector<double>> dNdY;
	resizeMatrixesForElementSize(H1, H2, dNdX, dNdY, Hbc1, Hbc2, c);

	for (int i = 0; i < grid.getNumberOfelements(); i++)
	{
		//cout << "Element id: " << i + 1 << endl;
		for (int j = 0; j < c.dNi.size(); j++)
		{

			inversejacobian = jacobian(i, j, grid, c, detJ);
			//	cout << "detj h: " << detJ << endl;;
			dNdxXdNdY(dNdX, dNdY, inversejacobian, c, j);
			calculatePieceOfHmatrix(dNdX[j], dNdY[j], H1, detJ, grid.elements[i].k);
			addMatrixH(H1, H2);

		}
		grid.elements[i].matrixH = H2;
		zeroMatrixH(H2);
	}
	cout << "inverse jacobian " << detJ << endl;
	zeroMatrixH(Hbc2);
	zeroMatrixH(Hbc1);
	for (int i = 0; i < grid.getNumberOfelements(); i++)
	{
		bool edgeElement = false;
		for (int j = 0; j < 4; j++)
		{
			if (grid.getNodeFromElement(i, j).BC) {
				edgeElement = true;
				break;
			}
		}
		if (edgeElement) {

			//cout << "det przy hbc: " << detJ<<endl;
			calculateHbc(Hbc1, Hbc2, c, grid.elements[i].alfa, grid, i);


		}
		else
			zeroMatrixH(Hbc2);


		grid.elements[i].matrixHBC = Hbc2;
		zeroMatrixH(Hbc2);

	}
	vector<double> p1, p2;

	makeVectorOfsize(p1, 4);
	makeVectorOfsize(p2, 4);
	for (int i = 0; i < grid.getNumberOfelements(); i++)
	{
		bool edgeElement = false;
		for (int j = 0; j < 4; j++)
		{
			if (grid.getNodeFromElement(i, j).BC) {
				edgeElement = true;
				break;
			}
		}
		if (edgeElement) {

			//cout << "det przy hbc: " << detJ<<endl;
			calculateVectorP(p1, p2, c, grid.elements[i].alfa, grid, i, grid.elements[i].ambient_temp);


		}
		else
			zeroVector(p2);


		grid.elements[i].vectorP = p2;
		zeroVector(p2);

	}



	vector<pair<double, double>> wspol;

	vector<pair<double, double>> points = makePoints(c, wspol);
	vector<vector<double>> C1;
	vector<vector<double>> C2;
	makeSquareMatrixOfSize(C1, 4);
	makeSquareMatrixOfSize(C2, 4);
	for (int i = 0; i < grid.getNumberOfelements(); i++)
	{
		//cout << "Element id: " << i + 1 << endl;
		for (int j = 0; j < c.dNi.size(); j++)
		{

			jacobian(i, j, grid, c, detJ);
			// cout << "detj c: " << detJ << endl;;
			auto shapeFun = calculatShapeFunctions(points[j]);

			calculatePeaceOfPieceOfC(C1, shapeFun, wspol[j], detJ, grid.elements[i].c, grid.elements[i].ro);

			addMatrixH(C1, C2);
			zeroMatrixH(C1);
		}
		grid.elements[i].matrixC = C2;
		zeroMatrixH(C2);
	}
	cout << "Matrix C: " << endl;
	grid.elements[0].printMatrixC();
	//for (auto x : grid.elements)
	//{
	//	x.printMatrixC();
	//	cout << "-----------------------------------------------------------" << endl;
	//}

	//cout << "detJ " << detJ<<endl;
	/*cout << "Matrix H: " << endl;
	for (auto x : grid.elements)
	{
		x.printMatrixH();
		cout << "-----------------------------------------------------------" << endl;
	}*/

	cout << "-----------------------------------------------------------" << endl;
	//cout << "Matrix Hbc:" << endl;
	//grid.elements[0].printMatrixHBC();
	/*for (auto x : grid.elements)
	{
		x.printMatrixHBC();
		cout << "-----------------------------------------------------------" << endl;
	}*/

	/*for (int k = 0; k < inversejacobian.size(); k++)
	{
		for (int l = 0; l < inversejacobian.size(); l++)
			cout << inversejacobian[k][l] << "     ";
		cout << endl;
	}*/
	/*cout << "--------------------------------\n";
	cout << "Print Vector p: \n";
	for (auto x : grid.elements)
	{
		x.printVectorP();
		cout << "-----------------------------------------------------------" << endl;
	}*/
	ThingsToCount thingsToCount(grid.nB * grid.nH);

	aggregateMatrixes(thingsToCount, grid);
	aggregateVectorP(thingsToCount, grid);
	thingsToCount.printAggregatedHmatrix();
	//thingsToCount.printAggregatedVectorP();
	thingsToCount.addHandC();
	thingsToCount.printAggregatedVectorP();

	/*thingsToCount.vectorPequation();
	solveAggregatedHmatrix(thingsToCount.aggregatedCmatrixPlusH, thingsToCount.solvedAggregatedHmatrix, thingsToCount.aggregatedVectorP, thingsToCount.equationSolution);




	thingsToCount.printEquationSolution();*/
	for (int i = 50; i <= 150; i += 50)    //here
	{ 
		thingsToCount.vectorPequation();
		solveAggregatedHmatrix(thingsToCount.aggregatedCmatrixPlusH, thingsToCount.solvedAggregatedHmatrix, thingsToCount.aggregatedVectorP, thingsToCount.equationSolution);


		//thingsToCount.printAggregatedVectorP();
		auto mx = max_element(std::begin(thingsToCount.equationSolution), std::end(thingsToCount.equationSolution));
		auto mn = min_element(std::begin(thingsToCount.equationSolution), std::end(thingsToCount.equationSolution));
		cout << "Time: " << i << "s min temp = " << *mn << " max temp = " << *mx << endl;
		//	thingsToCount.printEquationSolution();
	}
	//thingsToCount.printWholeSolvedAggregatedHmatrix();
	//thingsToCount.printAggregatedHmatrix();
	//	//Gauss(a, 1);
	//	//Gauss(b, 1);
	//	
	////	c.print();

	return 0;


}
