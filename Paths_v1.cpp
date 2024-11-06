#include<iostream>
#include<OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include <string>
#include <vector>
#include <algorithm>
#include<stdio.h>
#include "clipper.hpp"  
#include"clipper.cpp"
#include<queue>
#define use_int32
using namespace ClipperLib;
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
#define MAX 0.6
#define MIN 0.2
#define H 0.4  //默认的路径厚度（适用于Case2）
#define Width_path 0.8 //路径宽度，可以通过这个来调整填充率
#define Extrusion_width 0.8 //实际挤出的路径宽度
#define Width_material 1.75
#define pi 3.14
#define N 1	//首层与末层的层数
#define NP 4	//隔点上升

#define V_size 1.0 //用于支撑算法的体素大小
#define With_support 0.4 //支撑的路径宽度

using namespace std;
typedef unsigned char boolean;
struct Point {
	double x, y, z, t;
	bool b = false;
};

vector<Point>coord;  //缓存一层的轮廓线

vector<Point> buffer1;
vector<vector<Point>>buffer2;
vector<vector<vector<Point>>>paths;
vector<vector<vector<Point>>>P;

vector<vector<vector<Point>>>P_v;//支撑算法的体素点
vector<vector<vector<Point>>>Path_support;//支撑算法的路径

vector<vector<Point>>model; //slice  整个模型的，分层的轮廓线坐标，一维层数，二维此层的轮廓线点坐标
//vector<vector<Point>>model2; //slice  整个模型的，分层的轮廓线坐标，一维层数，二维此层的轮廓线点坐标
vector<Paths>model3;//从切片（model）转换类型得到的，做偏移之后画外轮廓。
vector<Point>shellpath;//存外壳上的点
vector<vector<Point> >shellpaths;//存某一层外壳上的路径
vector<vector<vector<Point> > >shell1;//存最里层轮廓上的路径      一维层数，二维每层不同的截面路径，三维路径的点坐标

vector<Paths>model2;//从切片（model）转换类型得到的，做偏移之后画外轮廓。
vector<vector<vector<Point> > >shell2;     //存某一层外壳上的路径集合  一维该层变化线数量，二维壳的层数路径，三维路径坐标


vector<vector<vector<vector<Point>>>> shell;//存所有外壳上的路径    一维层数，二维每层不同的截面路径，三维不同截面路径壳的层数，四维壳路径点的坐标


typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;//定义完整MyMesh类，用三角网格，自定的属性
MyMesh mesh;

double Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z, px, py, pz;//用于生成包围盒
int X, Y, Z;//包围盒内离散点的数量
int Scheme;

const string file_1 = "screw_big.stl";

// 读取文件的函数
void readfile(string file) {//文件读取到了mesh当中
	// 请求顶点法线 vertex normals
	mesh.request_vertex_normals();
	//如果不存在顶点法线，则报错 
	if (!mesh.has_vertex_normals())
	{
		cout << "错误：标准定点属性 “法线”不存在" << endl;
		return;
	}
	// 如果有顶点发现则读取文件 
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "无法读取文件:" << file << endl;
		return;
	}
	else cout << "成功读取文件:" << file << endl;
	cout << endl; // 为了ui显示好看一些
				  //如果不存在顶点法线，则计算出
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// 通过面法线计算顶点法线
		mesh.request_face_normals();
		// mesh计算出顶点法线
		mesh.update_normals();
		// 释放面法线
		mesh.release_face_normals();
	}
}

//截面函数
bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm) //pt截面上一点，pnorm截面法向
{
	const float ERR = 0.001;
	//参数包括 pt，pnorm，*pilist，pnum[]  具体函数原理 见 截面算法.docx
	int starte, ne, ne1, nf;
	MyMesh::Point vt1, vt2;
	//MyMesh::Face f1;
	MyMesh::HalfedgeHandle nhe;
	MyMesh::FaceHandle nhf;
	float d1, d2, sd1, sd2;
	bool* flag, suc;
	float dist, mind = 1.0e+8;

	sd1 = sd2 = -10000;
	int esize = mesh.n_halfedges();
	flag = new bool[esize];
	Point p;
	suc = false;


	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //遍历整个模型所有的边，有交点的把id记录在flag中
	{

		MyMesh::HalfedgeHandle hh = *it;
		int id = hh.idx();
		flag[id] = false;

		auto fromVertex = mesh.from_vertex_handle(hh);
		auto toVertex = mesh.to_vertex_handle(hh);
		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);
		//printf("$ %.3f %.3f $\n", vt1.data()[0],vt2.data()[0]);
		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);//d1到面的距离
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0) || d1 > 0 && d2 == 0 || d1 == 0 && d2 > 0))//线段与面相交
		{
			flag[id] = true;

			vt1.data()[0] = vt1.data()[0] - pt.data()[0];
			vt1.data()[1] = vt1.data()[1] - pt.data()[1];
			vt1.data()[2] = vt1.data()[2] - pt.data()[2];       // point date minus point date 
			dist = vt1.data()[0] * vt1.data()[0] + vt1.data()[1] * vt1.data()[1] + vt1.data()[2] * vt1.data()[2];

			if (dist < mind)
			{
				nhe = hh;	//最短边
				mind = dist;//最短距离
				ne = id;    //最短所在边的编号               //  printf("ne:  %d  \n", ne);
				suc = true;
			}
		}
	}

	if (!suc)
	{
		delete[]flag;
		return false; //没有交点，这里return false，跳出整个函数
	}

	starte = ne;//标记循环起始的边

	suc = false;

	nhf = mesh.face_handle(nhe);//最短边所在面

	while (!suc)
	{
		//printf("%%%%");	

		auto fromVertex = mesh.from_vertex_handle(nhe);
		auto toVertex = mesh.to_vertex_handle(nhe);

		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);

		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);
		//printf("$$$%lf %lf \n", d1, d2);
		if ((sd1 == d1) && (sd2 == d2))
		{
			flag[ne] = false;
		}
		sd1 = d1; sd2 = d2;
		//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
		//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
		//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];

		p = { (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[0] - vt1.data()[0]) + vt1.data()[0] ,(float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[1] - vt1.data()[1]) + vt1.data()[1] };
		coord.push_back(p);


		do {
			for (auto it = mesh.fh_begin(nhf); it != mesh.fh_end(nhf); ++it) //nhf最短边所在面,迭代这个面的边，只有3个
			{
				MyMesh::HalfedgeHandle halfnow = *it;

				const int ne1 = halfnow.idx();

				if (flag[ne1] == false || ne == ne1) continue;

				MyMesh::VertexHandle fromV = mesh.from_vertex_handle(halfnow);
				MyMesh::VertexHandle toV = mesh.to_vertex_handle(halfnow);

				vt1 = mesh.point(fromV);
				vt2 = mesh.point(toV);

				d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
				d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

				p = { (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[0] - vt1.data()[0]) + vt1.data()[0] ,(float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[1] - vt1.data()[1]) + vt1.data()[1] };
				coord.push_back(p);

				MyMesh::HalfedgeHandle halfnext = mesh.opposite_halfedge_handle(halfnow);//获取反向半边

				nhf = mesh.face_handle(halfnext);//返回这个边所在的面

				int ne2 = halfnext.idx();

				flag[ne1] = flag[ne2] = false;//ne1,ne2是对向的两个半边，存过其中一个的交点就都变为false

				if (nhf.idx() == -1)
				{
					starte = ne;
					flag[ne] = false;
					break;
				}
				ne = ne2;//以对边的那个面再找下一个与面相交的线

				//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
				//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
				//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
				//printf("##%lf %lf %lf\n", pilist[num].data()[0], pilist[num].data()[1], pilist[num].data()[2]);

				break;
			}
		} while (ne != starte);

		suc = true;

		for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //检索有没有第二个环
		{

			MyMesh::HalfedgeHandle hh = *it;

			int id = hh.idx();

			if (flag[id] == true)
			{
				ne = id;
				starte = ne;
				nhe = hh;
				nhf = mesh.face_handle(nhe);
				if (nhf.idx() == -1)
				{
					flag[ne] = false;
					continue;
				}
				//pilist[num].data()[0] = -10000;
				//pilist[num].data()[1] = -10000;
				//pilist[num].data()[2] = -10000;
				p = { -10000,-10000 };//两个环中间的间隔数据
				coord.push_back(p);

				suc = false;
				break;
			}
		}


	};

	delete[]flag;
	return true;
}


//通过遍历点获取模型的包围盒
void BoundingBox() {
	MyMesh::Point pt;
	int st = 0;
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		MyMesh::VertexHandle vh_i = *it;
		pt = mesh.point(vh_i);
		px = pt.data()[0];
		py = pt.data()[1];
		pz = pt.data()[2];
		if (st == 0) {
			Bmax_x = Bmin_x = px;
			Bmax_y = Bmin_y = py;
			Bmax_z = Bmin_z = pz;
			st++;
		}
		else {
			if (px > Bmax_x)Bmax_x = px; else if (px < Bmin_x)Bmin_x = px;
			if (py > Bmax_y)Bmax_y = py; else if (py < Bmin_y)Bmin_y = py;
			if (pz > Bmax_z)Bmax_z = pz; else if (pz < Bmin_z)Bmin_z = pz;
		}
	}
	Bmin_x -= Width_path / 2;
	Bmin_y -= Width_path / 2;
	//画图用
	/*Bmin_x -= Width_path / 2;
	Bmin_y -= Width_path / 2;
	Bmax_y += Width_path;
	Bmax_x += Width_path;*/

	X = (Bmax_x - Bmin_x) / Width_path;
	Y = (Bmax_y - Bmin_y) / Width_path;

	
}


double distance(Point a, Point b)	//两点间距离
{
	return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2) + pow((a.z - b.z), 2));
}

//示例一路径生成
void Case1_paths() {
	Point s;
	//生成离散点，距离为路径宽度
	Scheme = 1;
	s.x = Bmin_x; s.y = Bmin_y; s.z = 0; s.t = 0;
	Z = (Bmax_z - Bmin_z) / (MAX / 2 + MIN / 2) - 1;

	for (int i = 0; i < Z; i++) {
		if (i == 0) {//第一层
			for (int j = 0; j < X; j++) {
				s.x += Width_path;
				for (int m = 0; m < Y; m++) {
					s.y += Width_path;

					(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = Bmin_z + MIN / 2 : s.z = Bmin_z + MAX / 2;

					buffer1.push_back(s);
				}
				buffer2.push_back(buffer1);
				buffer1.clear();
				s.y = Bmin_y;
			}
			s.x = Bmin_x;
			P.push_back(buffer2);
			buffer2.clear();
		}
		else if (i == Z - 1) {//最后一层
			if ((i + 1) % 2 != 0) {//假如为奇数层
				for (int j = 0; j < X; j++) {
					s.x += Width_path;
					for (int m = 0; m < Y; m++) {
						s.y += Width_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN / 2 : s.z = P[i - 1][m][j].z + MAX / 2;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += Width_path;
					for (int m = 0; m < X; m++) {
						s.x += Width_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN / 2 : s.z = P[i - 1][m][j].z + MAX / 2;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
		else {//中间层
			if ((i + 1) % 2 != 0) {//假如为奇数层
				for (int j = 0; j < X; j++) {
					s.x += Width_path;
					for (int m = 0; m < Y; m++) {
						s.y += Width_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += Width_path;
					for (int m = 0; m < X; m++) {
						s.x += Width_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
	}

}

//示例二路径生成
void Case2_paths() {
	Point s;
	Scheme = 2;
	s.x = Bmin_x; s.y = Bmin_y; s.z = 0; s.t = 0;
	Z = (Bmax_z - Bmin_z) / (MAX / 2 + MIN / 2) - 1;
	for (int i = 0; i < Z; i++) {
		if (i < N) {//首层
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += Width_path;//假如为奇数层
					for (int m = 0; m < Y; m++) {
						s.y += Width_path;
						
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = Bmin_z + MIN : s.z = Bmin_z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += Width_path;//假如为偶数层
					for (int m = 0; m < X; m++) {
						s.x += Width_path;
					
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
		else if (i >= Z - N) {//末层//假如为奇数层
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += Width_path;//假如为奇数层
					for (int m = 0; m < Y; m++) {
						s.y += Width_path;
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MAX : s.z = P[i - 1][m][j].z + MIN;//不用考虑奇偶层，应该与第一层互补
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += Width_path;
					for (int m = 0; m < X; m++) {
						s.x += Width_path;
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MAX : s.z = P[i - 1][m][j].z + MIN;//不用考虑奇偶层，应该与第一层互补
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
		else {//中间层//假如为奇数层
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += Width_path;
					for (int m = 0; m < Y; m++) {
						s.y += Width_path;
						s.z = P[i - 1][m][j].z + H;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += Width_path;
					for (int m = 0; m < X; m++) {
						s.x += Width_path;
						s.z = P[i - 1][m][j].z + H;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}

	}
}

//判断点是否在多边形内，在true，否false，使用的是引射线法
bool InOrOutPolygon(Point a, vector<Point> polypoint) {
	double x0 = a.x;
	double y0 = a.y;
	int crossings = 0;
	int n = polypoint.size();
	for (int i = 0; i < n; i++)
	{
		// 点在两个x之间 且以点垂直y轴向上做射线
		double slope = (polypoint[(i + 1 + n) % n].y - polypoint[i].y) / (polypoint[(i + 1 + n) % n].x - polypoint[i].x);
		boolean cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		boolean cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		boolean above = (y0 < slope* (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	if (crossings % 2 != 0&&crossings!=0) {
		return true;
	}
	else {
		return false;
	}
}

//判断生成的离散点是否在模型内
void findIntersect() {
	Path  way;
	Paths layer;//每层路径集合
	IntPoint p1;
	Point p;
	for (int i = 0; i < P.size(); i++) {
		printf("求交算法：%.2lf%%\r", i * 100.0 / P.size());
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 0; m < P[i][j].size(); m++) {
				MyMesh::Normal vf(0, 0, 1);//截面法线
				MyMesh::Point pt;//截面上一点
				pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = P[i][j][m].z;
				IntersectPlane(pt, vf);//生成一层轮廓线，存入coord

				for (int k = 0; k < coord.size(); k++)
				{
					if (coord[k].x != -10000)
					{
						p1.X = coord[k].x * 100;
						p1.Y = coord[k].y * 100;
						way.push_back(p1);
					}
					else
					{
						layer.push_back(way);
						way.clear();
					}
				}
				layer.push_back(way);
				way.clear();
				coord.clear();

				for (int k1 = 0; k1 < layer.size(); k1++)
				{
					for (int k = 0; k < layer[k1].size(); k++)
					{
						p.x = double(layer[k1][k].X) / 100;
						p.y = double(layer[k1][k].Y) / 100;

						shellpath.push_back(p);
					}
					shellpaths.push_back(shellpath);
					shellpath.clear();
				}
				layer.clear();

				for (int k = 0; k < shellpaths.size(); k++) {
					if (InOrOutPolygon(P[i][j][m], shellpaths[k])) {
						P[i][j][m].b = true;
					}
				}
				shellpaths.clear();
				
			}
		}
		
	}
}

void Contour() {
	double z = Bmin_z;
	for (int i = 0; i < Z; i++) {
		z += (MIN + MAX) / 2;
		MyMesh::Normal vf(0, 0, 1);//截面法向
		MyMesh::Point pt; //截面上一点
		pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = z;
		IntersectPlane(pt, vf);//生成一层轮廓线，存入coord
		model.push_back(coord);
		coord.clear();
	}
	Path  way;
	Paths layer;//每层路径集合
	IntPoint p1;
	ClipperOffset co;
	Paths solution;
	Point p;
	for (int i = 0; i < model.size(); i++)//遍历model，第i层
	{
	

		for (int j = 0; j < model[i].size(); j++)//遍历model，第i层中的第j个点
		{
			if (model[i][j].x != -10000)
			{
				p1.X = model[i][j].x * 100;
				p1.Y = model[i][j].y * 100;
				way.push_back(p1);
			}
			else
			{
				layer.push_back(way);
				way.clear();
			}
		}
		layer.push_back(way);
		way.clear();

		model2.push_back(layer);
		layer.clear();
//*********************将最内圈的点赋给shell1***********************************
		for (int j = 0; j < model2[i].size(); j++)
		{
			for (int k = 0; k < model2[i][j].size(); k++)
			{
				p.x = double(model2[i][j][k].X) / 100;
				p.y = double(model2[i][j][k].Y) / 100;

				shellpath.push_back(p);
			}
			shellpaths.push_back(shellpath);
			shellpath.clear();
		}
		shell1.push_back(shellpaths);
		shellpaths.clear();



		//******************************做偏移***********************************************
		printf("%.2lf%%\r", i * 100.0 / model.size());
		for (int j = 0; j < model2[i].size(); j++)//此时遍历的是model2，某层里的第j条path做偏移
		{
			double area = abs(Area(model2[i][j]) / 10000);
			co.Clear();
			co.AddPath(model2[i][j], jtRound, etClosedPolygon);   //设置准备偏移的路径
			int times = 1;	//偏移次数，即外壳打印的层数
			int sh = Extrusion_width * 100;     //偏移厚度

			for (int a = 0; a < times; a++)					 //进行偏移
			{

				co.Execute(solution, sh);							//每次偏移的厚度

				for (int n = 0; n < solution[0].size(); n++)
				{
					p.x = double(solution[0][n].X) / 100;
					p.y = double(solution[0][n].Y) / 100;

					shellpath.push_back(p);
				}

				//p.x = double(solution[0][0].X) / 100;
				//p.y = double(solution[0][0].Y) / 100;

				shellpath.push_back(p);

				co.Clear();
				co.AddPath(solution[0], jtRound, etClosedPolygon);  //准备下一次偏移
				solution.clear();

				shellpaths.push_back(shellpath);
				shellpath.clear();
			}

			shell2.push_back(shellpaths);
			shellpaths.clear();

		}
		//solution.clear();
		shell.push_back(shell2);
		shell2.clear();
	}

}

void Optimization() {
	Point s;
	//double k1,k2;	//斜率
	double g;//插值间隔
	for (int i = 0; i < P.size(); i++) {//中间插入过渡插值
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {
				if (P[i][j][m].z != P[i][j][m - 1].z) {
					//k1 = (P[i][j][m].z - P[i][j][m - 1].z) / Width_path;
					//k2 = (P[i][j][m].t - P[i][j][m - 1].t) / Width_path;
					Scheme == 1 ? g = 0.025 : g = 0.05;
					for (int c = 1; c < 8; c++) {//
						P[i][j][m - 1].x == P[i][j][m].x ? s.y = P[i][j][m - 1].y + Width_path / 8 : s.x = P[i][j][m - 1].x + Width_path / 8;
						P[i][j][m - 1].x == P[i][j][m].x ? s.x = P[i][j][m - 1].x : s.y = P[i][j][m - 1].y;
						P[i][j][m].z > P[i][j][m - 1].z ? s.z = P[i][j][m - 1].z + g : s.z = P[i][j][m - 1].z - g;
						//s.z = P[i][j][m - 1].z + g * k1;
						//s.t = P[i][j][m - 1].t + 0.1 * k2;
						(P[i][j][m - 1].b == true && P[i][j][m].b == true) ? s.b = true : s.b = false;
						P[i][j].insert(P[i][j].begin() + m, s);
						m += 1;
					}
				}

			}
		}
	}

}

void  Case1() {
	Case1_paths();
	findIntersect();  //前三个函数切片，产生model
	Optimization();

	double k = (MAX / 2 - MIN / 2) / Width_path;//存储斜率
	double z = (MIN + MAX) / 2;
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {

				//首层
				if (i == 0) {
					P[i][j][m].t = 2 * (Extrusion_width) * (P[i][j][m].z + P[i][j][m - 1].z) / (pow(Width_material, 2) * pi);
				}
				else if (i == 1) {//第二层
					if (P[i][j][m].z != P[i][j][m - 1].z) {
						P[i][j][m].t = 2 * (Extrusion_width) * 2 *(P[i][j][m].z - 0.5 + MIN / 2 + P[i][j][m - 1].z - 0.5 + MIN / 2) / (pow(Width_material, 2) * pi);
					}
					/*if (P[i][j][m].z > P[i][j][m - 1].z) {
						
						P[i][j][m].t = 2 * (Extrusion_width) * (P[i][j][m].z - 0.5 + MIN / 2 + P[i][j][m - 1].z - 0.5 + MIN / 2) / (pow(Width_material, 2) * pi);

					}
					else if (P[i][j][m].z < P[i][j][m - 1].z) {
						
							P[i][j][m].t = 2 * (Extrusion_width) * ((MAX - 2 * ((m % 8) * (Width_path / 8) * k)) + (MAX - 2 * (((m - 1) % 8) * (Width_path / 8) * k))) / (pow(Width_material, 2) * pi);
							 
					}*/
					else {
						if (P[i][j][m].z == 0.7)
							P[i][j][m].t = 4 * (Extrusion_width)*MAX / (pow(Width_material, 2) * pi);
						else
							P[i][j][m].t = 4 * (Extrusion_width)*MIN / (pow(Width_material, 2) * pi);
					}
				}//末层(奇数)
				else if (i == Z - 1 && (i + 1) % 2 != 0) {
					P[i][j][m].t = P[0][j][m].t;
				}//末层（偶数）
				else if (i == Z - 1 && (i + 1) % 2 == 0) {
					P[i][j][m].t = P[1][j][m].t / 2;
				}//中间层
				else {
					if ((i + 1) % 2 != 0) {
						P[i][j][m].t = P[0][j][m].t * 2;
					}
					else {
						P[i][j][m].t = P[1][j][m].t;
					}

				}

			}
		}
	}
}

void  Case2() {
	Case2_paths();
	findIntersect();  //前三个函数切片，产生model
	Optimization();

	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {
				//首层
				if (i < N) {
					if (i == 0) {
						P[i][j][m].t = 2 * (Extrusion_width) * (P[i][j][m].z + P[i][j][m - 1].z) / (pow(Width_material, 2) * pi);
					}
					else {
						P[i][j][m].t = 2 * (Extrusion_width) * (P[i][j][m].z/2 + P[i][j][m - 1].z/2) / (pow(Width_material, 2) * pi);
					}
				}
				//末层(奇数)
				else if (i >= Z - N && (i + 1) % 2 != 0) {
					P[i][j][m].t = 2 * (Extrusion_width) * (MAX + MIN - P[0][j][m - 1].z + MAX + MIN - P[0][j][m].z) / (pow(Width_material, 2) * pi);
				} // 末层(偶数)
				else if (i >= Z - N && (i + 1) % 2 == 0) {
					P[i][j][m].t = 2 * (Extrusion_width) * (MAX + MIN - (P[1][j][m - 1].z - H) + MAX + MIN - (P[1][j][m].z - H)) / (pow(Width_material, 2) * pi);

				}//中间层
				else {
					P[i][j][m].t = 4 * (Extrusion_width) * H / (pow(Width_material, 2) * pi);
				}
			}
		}
	}
}

void Paths_Optimization() {
	Point s;
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 0; m < P[i][j].size(); m++) {
				if (P[i][j][m].b == true) {
					buffer1.push_back(P[i][j][m]);
				}
				else if (P[i][j][m].b == false && buffer1.size() != 0) {
					buffer2.push_back(buffer1);
					buffer1.clear();
				}
			}
			if (buffer1.size() != 0) {
				buffer2.push_back(buffer1);
				buffer1.clear();
			}
		}
		paths.push_back(buffer2);
		buffer2.clear();
	}

	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			if (i % 2 == 0) {
				paths[i][j][0].y -= Width_path / 2;
				paths[i][j][paths[i][j].size() - 1].y += Width_path / 2;
			}
			else {
				paths[i][j][0].x -= Width_path / 2;
				paths[i][j][paths[i][j].size() - 1].x += Width_path / 2;
			}
		}
	}

	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size(); j++) {
			if (j % 2 == 0) {
				reverse(paths[i][j].begin(), paths[i][j].end());
			}
		}
	}




	double arcs, s_begin, s_end;
	for (int i = 0; i < paths.size()-10; i++) {
		for (int j = 0; j < paths[i].size() - 1; j++) {
			arcs = distance(paths[i][j][paths[i][j].size() - 1], paths[i][j + 1][0]);
			for (int m = j + 1; m < paths[i].size(); m++) {
				s_begin = distance(paths[i][j][paths[i][j].size() - 1], paths[i][m][0]);
				s_end = distance(paths[i][j][paths[i][j].size() - 1], paths[i][m][paths[i][m].size() - 1]);
				if (s_begin < arcs && s_begin < s_end) {
					arcs = s_begin;
					buffer1 = paths[i][m];
					paths[i][m] = paths[i][j + 1];
					paths[i][j + 1] = buffer1;
					buffer1.clear();
				}
				else if (s_end < arcs && s_end < s_begin) {
					arcs = s_end;
					reverse(paths[i][m].begin(), paths[i][m].end());
					buffer1 = paths[i][m];
					paths[i][m] = paths[i][j + 1];
					paths[i][j + 1] = buffer1;
					buffer1.clear();
				}
			}
		}
	}

	

}

void ContourGcodePrint() {


	FILE* fp;
	errno_t err;     //判断此文件流是否存在 存在返回1
	err = fopen_s(&fp, "ls2.gcode", "a"); //若return 1 , 则将指向这个文件的文件流给
	double t = (4 * ((MIN + MAX) / 2) * Extrusion_width) / (pow(Width_material, 2) * pi);
	double E = 0;
	double r;//回抽
	int L = 50;//偏移量
	fprintf(fp, ";FLAVOR:Marlin\n");
	fprintf(fp, ";Generated with Cura_SteamEngine 4.10.0\n");
	fprintf(fp, "M140 S50\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M190 S50\n");
	fprintf(fp, "M104 S210\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S210\n");
	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M201 X500.00 Y500.00 Z100.00 E5000.00 ;Setup machine max acceleration\n");
	fprintf(fp, "M203 X500.00 Y500.00 Z10.00 E50.00 ;Setup machine max feedrate\n");
	fprintf(fp, "M204 P500.00 R1000.00 T500.00 ;Setup Print/Retract/Travel acceleration\n");
	fprintf(fp, "M205 X8.00 Y8.00 Z0.40 E5.00 ;Setup Jerk\n");
	fprintf(fp, "M220 S100 ;Reset Feedrate\n");
	fprintf(fp, "M221 S100 ;Reset Flowrate\n");

	fprintf(fp, "G28 ;Home\n");

	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F2700 E-5\n");
	fprintf(fp, "M107\n");

	double z = Bmin_z;

	for (int i = 0; i < shell1.size(); i++) {
		fprintf(fp, ";LAYER:%d\n", i);
		fprintf(fp, "M73 P%.f\n", float((100 * i) / paths.size()));
		fprintf(fp, ";TYPE:OUTLINE\n");
		z += (MIN + MAX) / 2;
		fprintf(fp, "G0 F4800 X%f Y%f Z%f\n", shell1[i][0][0].x + L, shell1[i][0][0].y + L, z);
		if (i > 0) {
			fprintf(fp, "G1 F2400 E%f\n", E += 1);
		}
		for (int j = 0; j < shell1[i].size(); j++) {
			fprintf(fp, "G0 X%f Y%f Z%f\n", shell1[i][j][0].x + L, shell1[i][j][0].y + L, z);
			 
			fprintf(fp, "G1 F1200 X%f Y%f Z%f E%f\n", shell1[i][j][1].x + L, shell1[i][j][1].y + L,z, E += distance(shell1[i][j][0], shell1[i][j][1]) * t);
			for (int k = 2; k < shell1[i][j].size(); k++)

			{
				fprintf(fp, "G1 X%.3f Y%.3f Z%f E%.5f\n", shell1[i][j][k].x + L, shell1[i][j][k].y + L,z, E += distance(shell1[i][j][k - 1], shell1[i][j][k]) * t);

			}
			for (int k = 0; k < shell[i][j].size(); k++)
			{
				fprintf(fp, "G0 F4800 X%.3f Y%.3f Z%.1f\n", shell[i][j][k][0].x + L, shell[i][j][k][0].y + L, z);

				fprintf(fp, "G1 F1200 X%.3f Y%.3f Z%f E%.5f\n", shell[i][j][k][1].x + L, shell[i][j][k][1].y + L,z, E += distance(shell[i][j][k][0], shell[i][j][k][1]) * t);

				for (int m = 2; m < shell[i][j][k].size(); m++)

				{
					fprintf(fp, "G1 X%.3f Y%.3f Z%f E%.5f\n", shell[i][j][k][m].x + L, shell[i][j][k][m].y + L,z, E += distance(shell[i][j][k][m - 1], shell[i][j][k][m]) * t);

				}
				fprintf(fp, "G1 F1200 X%.3f Y%.3f Z%.1f\n", shell[i][j][k][0].x + L, shell[i][j][k][0].y + L, z, E += distance(shell[i][j][k][shell[i][j][k].size() - 1], shell[i][j][k][0]) * t);
			}
		}
		
		fprintf(fp, ";TYPE:FILL\n");
		if (i <paths.size() && paths.size() != 0) {
			for (int j = 0; j < paths[i].size(); j++) {
				r = E - 1;
				fprintf(fp, "G1 F2400 E%f\n", r);//回抽一定距离，避免拉丝
				fprintf(fp, "G0 X%f Y%f Z%f\n", paths[i][j][0].x + L, paths[i][j][0].y + L, paths[i][j][0].z);
				fprintf(fp, "G1 F2400 E%f\n", E);
				for (int m = 1; m < paths[i][j].size(); m++) {

					fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", paths[i][j][m].x + L, paths[i][j][m].y + L, paths[i][j][m].z, E += sqrt(pow(abs(paths[i][j][m].x - paths[i][j][m - 1].x), 2) + pow(abs(paths[i][j][m].y - paths[i][j][m - 1].y), 2)) * paths[i][j][m].t);

				}
			}
		}

		fprintf(fp, ";TYPE:SUPPORT\n");
		if (i < Path_support.size() && Path_support[i].size() != 0) {
			for (int j = 0; j < Path_support[i].size(); j++) {
				r = E - 1;
				fprintf(fp, "G1 F2400 E%f\n", r);//回抽一定距离，避免拉丝

				fprintf(fp, "G0 F4800 X%f Y%f\n", Path_support[i][j][0].x + L, Path_support[i][j][0].y + L);
				fprintf(fp, "G1 F2400 E%f\n", E);

				for (int m = 1; m < Path_support[i][j].size(); m++) {

					fprintf(fp, "G1 F1200 X%f Y%f E%f\n", Path_support[i][j][m].x + L, Path_support[i][j][m].y + L, E += distance(Path_support[i][j][m - 1], Path_support[i][j][m]) * 0.04756);
				}
			}

		}

		fprintf(fp, "G0 F3000 Z%f E%f\n", paths[i][paths[i].size() - 1][paths[i][paths[i].size() - 1].size() - 1].z + MAX, E += -1);

	}
	fprintf(fp, "M140 S0\n");
	fprintf(fp, "M107\n");
	fprintf(fp, "G91\n");
	fprintf(fp, "G1 E-2 F2700\n");
	fprintf(fp, "G1 E-2 Z0.2 F2400 ;Retract and raise Z\n");
	fprintf(fp, "G1 X5 Y5 F3000 ;Wipe out\n");
	fprintf(fp, "G1 Z10 ;Raise Z more\n");
	fprintf(fp, "G90 ;Absolute positioning\n");

	fprintf(fp, "G1 X0 Y300 ;Present print\n");
	fprintf(fp, "M106 S0 ;Turn-off fan\n");
	fprintf(fp, "M104 S0 ;Turn-off hotend\n");
	fprintf(fp, "M140 S0 ;Turn-off bed\n");

	fprintf(fp, "M84 X Y E ;Disable all steppers but Z\n");

	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fclose(fp);
}

//void Support() {
//
//	Point p;
//	//创建包围盒体素块
//	for (double i = Bmin_z; i < Bmax_z; i += H) {
//		for (double j = Bmin_y + V_size / 2; j < Bmax_y + V_size / 2; j += V_size) {
//			for (double m = Bmin_x + V_size/2; m < Bmax_x + V_size/2; m += V_size) {
//				p.x = m; p.y = j; p.z = i;
//				p.t = 2 * With_support * V_size / (pow(Width_material, 2) * pi);
//				p.b = true;
//				buffer1.push_back(p);
//			}
//			buffer2.push_back(buffer1);
//			buffer1.clear();
//		}
//		P_v.push_back(buffer2);
//		buffer2.clear();
//	}
//
//	//判断轮廓外的点
//	Path  way;
//	Paths layer;//每层路径集合
//	IntPoint p1;
//	for (int i = 0; i < P_v.size(); i++) {
//		printf("求交算法：%.2lf%%\r", i * 100.0 / P_v.size());
//		for (int j = 0; j < P_v[i].size(); j++) {
//			for (int m = 0; m < P_v[i][j].size(); m++) {
//				MyMesh::Normal vf(0, 0, 1);//截面法线
//				MyMesh::Point pt;//截面上一点
//				pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = P_v[i][j][m].z;
//				IntersectPlane(pt, vf);//生成一层轮廓线，存入coord
//
//				for (int k = 0; k < coord.size(); k++)
//				{
//					if (coord[k].x != -10000)
//					{
//						p1.X = coord[k].x * 100;
//						p1.Y = coord[k].y * 100;
//						way.push_back(p1);
//					}
//					else
//					{
//						layer.push_back(way);
//						way.clear();
//					}
//				}
//				layer.push_back(way);
//				way.clear();
//				coord.clear();
//
//				for (int k1 = 0; k1 < layer.size(); k1++)
//				{
//					for (int k = 0; k < layer[k1].size(); k++)
//					{
//						p.x = double(layer[k1][k].X) / 100;
//						p.y = double(layer[k1][k].Y) / 100;
//
//						shellpath.push_back(p);
//					}
//					shellpaths.push_back(shellpath);
//					shellpath.clear();
//				}
//				layer.clear();
//
//				for (int k = 0; k < shellpaths.size(); k++) {
//					if (InOrOutPolygon(P_v[i][j][m], shellpaths[k])) {
//						P_v[i][j][m].b = false;
//					}
//				}
//				shellpaths.clear();
//
//			}
//		}
//
//	}
//	//生成支撑路径
//	Point s;
//	for (int i = 0; i < P_v.size(); i++) {
//		for (int j = 0; j < P_v[i].size(); j++) {
//			for (int m = 0; m < P_v[i][j].size(); m++) {
//				if (P_v[i][j][m].b == true) {
//					buffer1.push_back(P_v[i][j][m]);
//				}
//				else if (P_v[i][j][m].b == false && buffer1.size() != 0) {
//					buffer2.push_back(buffer1);
//					buffer1.clear();
//				}
//			}
//			if (buffer1.size() != 0) {
//				buffer2.push_back(buffer1);
//				buffer1.clear();
//			}
//		}
//		Path_support.push_back(buffer2);
//		buffer2.clear();
//	}
//
//	double arcs, s_begin, s_end;
//	for (int i = 0; i < Path_support.size(); i++) {
//		for (int j = 0; j < Path_support[i].size() - 1; j++) {
//			arcs = distance(Path_support[i][j][Path_support[i][j].size() - 1], Path_support[i][j + 1][0]);
//			for (int m = j + 1; m < Path_support[i].size(); m++) {
//				s_begin = distance(Path_support[i][j][Path_support[i][j].size() - 1], Path_support[i][m][0]);
//				s_end = distance(Path_support[i][j][Path_support[i][j].size() - 1], Path_support[i][m][Path_support[i][m].size() - 1]);
//				if (s_begin < arcs && s_begin < s_end) {
//					arcs = s_begin;
//					buffer1 = Path_support[i][m];
//					Path_support[i][m] = Path_support[i][j + 1];
//					Path_support[i][j + 1] = buffer1;
//					buffer1.clear();
//				}
//				else if (s_end < arcs && s_end < s_begin) {
//					arcs = s_end;
//					reverse(Path_support[i][m].begin(), Path_support[i][m].end());
//					buffer1 = Path_support[i][m];
//					Path_support[i][m] = Path_support[i][j + 1];
//					Path_support[i][j + 1] = buffer1;
//					buffer1.clear();
//				}
//			}
//		}
//	}
//
//
//}

void Support() {
	Point s;
	int m;
	int p2;
	int p3;

	Point p;
		//创建包围盒体素块
		for (double i = Bmin_z; i < Bmax_z; i += H) {
			for (double j = Bmin_y + V_size / 2; j < Bmax_y + V_size / 2; j += V_size) {
				for (double m = Bmin_x + V_size/2; m < Bmax_x + V_size/2; m += V_size) {
					p.x = m; p.y = j; p.z = i;
					p.t = 2 * With_support * V_size / (pow(Width_material, 2) * pi);
					p.b = false;
					buffer1.push_back(p);
				}
				buffer2.push_back(buffer1);
				buffer1.clear();
			}
			P_v.push_back(buffer2);
			buffer2.clear();
		}
	
		//判断轮廓外的点
		Path  way;
		Paths layer;//每层路径集合
		IntPoint p1;
		for (int i = 0; i < P_v.size(); i++) {
			printf("求交算法：%.2lf%%\r", i * 100.0 / P_v.size());
			for (int j = 0; j < P_v[i].size(); j++) {
				for (int m = 0; m < P_v[i][j].size(); m++) {
					MyMesh::Normal vf(0, 0, 1);//截面法线
					MyMesh::Point pt;//截面上一点
					pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = P_v[i][j][m].z;
					IntersectPlane(pt, vf);//生成一层轮廓线，存入coord
	
					for (int k = 0; k < coord.size(); k++)
					{
						if (coord[k].x != -10000)
						{
							p1.X = coord[k].x * 100;
							p1.Y = coord[k].y * 100;
							way.push_back(p1);
						}
						else
						{
							layer.push_back(way);
							way.clear();
						}
					}
					layer.push_back(way);
					way.clear();
					coord.clear();
	
					for (int k1 = 0; k1 < layer.size(); k1++)
					{
						for (int k = 0; k < layer[k1].size(); k++)
						{
							p.x = double(layer[k1][k].X) / 100;
							p.y = double(layer[k1][k].Y) / 100;
	
							shellpath.push_back(p);
						}
						shellpaths.push_back(shellpath);
						shellpath.clear();
					}
					layer.clear();
	
					for (int k = 0; k < shellpaths.size(); k++) {
						if (InOrOutPolygon(P_v[i][j][m], shellpaths[k])) {
							P_v[i][j][m].b = true;
						}
					}
					shellpaths.clear();
	
				}
			}
	
		}

	int c = 2; //c表示者周围的几个单元是否存在真值，属于一个范围区间，1代表着周围的八个个，2就更大。主要防止支撑与最外的轮廓相交，从而影响表面精度。
	for (int i = 0; i < P_v.size(); i++) {
		for (int j = c; j < P_v[i].size() - c; j++) {
			j % 2 == 0 ? m = c : m = P_v[i][j].size() - 1 - c;
			while (m > c-1 && m < P_v[i][j].size() - c) {

				if (P_v[i][j][m].b == false) {
					p2 = 0;
					p3 = i;
					while (p3 < P_v.size()) {
						if (P_v[p3][j][m].b == true) {
							p2 = 1;
							break;
						}
						p3++;
					}
					if (i > 0) {
						if (P_v[i][j][m + c].b == true || P_v[i][j][m - c].b == true || P_v[i][j + c][m].b == true || P_v[i][j - c][m].b == true || P_v[i][j + c][m + c].b == true || P_v[i][j + c][m - c].b == true
							|| P_v[i][j - c][m + c].b == true || P_v[i][j - c][m - c].b == true || P_v[i - 1][j][m].b == true) {
							p2 = 0;
						}
					}

					if (p2 == 1) {
						s.x = P_v[i][j][m].x; s.y = P_v[i][j][m].y;
						buffer1.push_back(s);
					}
					else if (p2 == 0 && buffer1.size() != 0) {
						buffer2.push_back(buffer1);
						buffer1.clear();
					}

				}
				else if (P_v[i][j][m].b == true && buffer1.size() != 0) {
					buffer2.push_back(buffer1);
					buffer1.clear();
				}

				j % 2 == 0 ? m++ : m--;
			}
			if (buffer1.size() != 0) {
				buffer2.push_back(buffer1);
				buffer1.clear();
			}
		}
		Path_support.push_back(buffer2);
		buffer2.clear();
	}

	double arcs, s_begin, s_end;
	for (int i = 0; i < Path_support.size(); i++) {
		for (int j = 1; j < Path_support[i].size() ; j++) {
			
			arcs = distance(Path_support[i][j - 1][Path_support[i][j - 1].size() - 1], Path_support[i][j][0]);
			for (int m = j; m < Path_support[i].size(); m++) {
				s_begin = distance(Path_support[i][j-1][Path_support[i][j - 1].size() - 1], Path_support[i][m][0]);
				s_end = distance(Path_support[i][j - 1][Path_support[i][j - 1].size() - 1], Path_support[i][m][Path_support[i][m].size() - 1]);
				if (s_begin < arcs && s_begin < s_end) {
					arcs = s_begin;
					buffer1 = Path_support[i][m];
					Path_support[i][m] = Path_support[i][j];
					Path_support[i][j] = buffer1;
					buffer1.clear();
				}
				else if (s_end < arcs && s_end < s_begin) {
					arcs = s_end;
					reverse(Path_support[i][m].begin(), Path_support[i][m].end());
					buffer1 = Path_support[i][m];
					Path_support[i][m] = Path_support[i][j];
					Path_support[i][j] = buffer1;
					buffer1.clear();
				}
			}
		}
	}




}


void GcodePrint() {
	FILE* fp;
	errno_t err;     //判断此文件流是否存在 存在返回1
	err = fopen_s(&fp, "ls2.gcode", "a"); //若return 1 , 则将指向这个文件的文件流给
	double t = (4 * ((MIN + MAX) / 2) * Width_path) / (pow(Width_material, 2) * pi);
	double E = 0;
	double r;//回抽
	int L = 50;//偏移量
	fprintf(fp, ";FLAVOR:Marlin\n");
	fprintf(fp, ";Generated with Cura_SteamEngine 4.10.0\n");
	fprintf(fp, "M140 S50\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M190 S50\n");
	fprintf(fp, "M104 S210\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S210\n");
	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M201 X500.00 Y500.00 Z100.00 E5000.00 ;Setup machine max acceleration\n");
	fprintf(fp, "M203 X500.00 Y500.00 Z10.00 E50.00 ;Setup machine max feedrate\n");
	fprintf(fp, "M204 P500.00 R1000.00 T500.00 ;Setup Print/Retract/Travel acceleration\n");
	fprintf(fp, "M205 X8.00 Y8.00 Z0.40 E5.00 ;Setup Jerk\n");
	fprintf(fp, "M220 S100 ;Reset Feedrate\n");
	fprintf(fp, "M221 S100 ;Reset Flowrate\n");

	fprintf(fp, "G28 ;Home\n");

	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F2700 E-5\n");
	fprintf(fp, "M107\n");

	double z = Bmin_z;

	for (int i = 0; i < paths.size(); i++) {
		fprintf(fp, ";LAYER:%d\n", i);
	
		fprintf(fp, "M73 P%.f\n", float((100 * i) / paths.size()));
		fprintf(fp, "G0 X%f Y%f Z%f\n", paths[i][0][0].x + L, paths[i][0][0].y + L, paths[i][0][0].z);
		fprintf(fp, "G0 E%f\n", E += 0.2);

		fprintf(fp, ";TYPE:FILL\n");
		for (int j = 0; j < paths[i].size(); j++) {
			
			fprintf(fp, "G0 X%f Y%f Z%f\n", paths[i][j][0].x + L, paths[i][j][0].y + L, paths[i][j][0].z);
			if (i > 0) {
				
			}
			for (int m = 1; m < paths[i][j].size(); m++) {

				fprintf(fp, "G1 F1200 X%f Y%f Z%f E%f\n", paths[i][j][m].x + L, paths[i][j][m].y + L, paths[i][j][m].z, E += sqrt(pow(abs(paths[i][j][m].x - paths[i][j][m - 1].x), 2) + pow(abs(paths[i][j][m].y - paths[i][j][m - 1].y), 2)) * paths[i][j][m].t);

			}

		}

		//support
		/*if (i > 0) {
			for (int j = 0; j < Path_support[i].size(); j++) {

				fprintf(fp, "G0 X%f Y%f Z%f\n", Path_support[i][j][0].x + L, Path_support[i][j][0].y + L, Path_support[i][j][0].z);
				if (i > 0) {

				}
				for (int m = 1; m < Path_support[i][j].size(); m++) {

					fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", Path_support[i][j][m].x + L, Path_support[i][j][m].y + L, Path_support[i][j][m].z, E += distance(Path_support[i][j][m - 1], Path_support[i][j][m]) * Path_support[i][j][m].t);

				}

			}
		}*/

		if (i < Path_support.size() && Path_support[i].size() != 0) {
			for (int j = 0; j < Path_support[i].size(); j++) {
				r = E - 0.2;
				fprintf(fp, "G1 F300 E%f\n", r);//回抽一定距离，避免拉丝

				fprintf(fp, "G0 F4800 X%f Y%f\n", Path_support[i][j][0].x + L, Path_support[i][j][0].y + L);
				fprintf(fp, "G1 F300 E%f\n", E);

				for (int m = 1; m < Path_support[i][j].size(); m++) {

					fprintf(fp, "G1 F1200 X%f Y%f E%f\n", Path_support[i][j][m].x + L, Path_support[i][j][m].y + L, E += distance(Path_support[i][j][m - 1], Path_support[i][j][m]) * 0.04756);
				}
			}

		}

		fprintf(fp, "G0 F3000 Z%f E%f\n", paths[i][paths[i].size() - 1][paths[i][paths[i].size() - 1].size() - 1].z + MAX, E -= 0.2);


	}
	fprintf(fp, "M140 S0\n");
	fprintf(fp, "M107\n");
	fprintf(fp, "G91\n");
	fprintf(fp, "G1 E-2 F2700\n");
	fprintf(fp, "G1 E-2 Z0.2 F2400 ;Retract and raise Z\n");
	fprintf(fp, "G1 X5 Y5 F3000 ;Wipe out\n");
	fprintf(fp, "G1 Z10 ;Raise Z more\n");
	fprintf(fp, "G90 ;Absolute positioning\n");

	fprintf(fp, "G1 X0 Y300 ;Present print\n");
	fprintf(fp, "M106 S0 ;Turn-off fan\n");
	fprintf(fp, "M104 S0 ;Turn-off hotend\n");
	fprintf(fp, "M140 S0 ;Turn-off bed\n");

	fprintf(fp, "M84 X Y E ;Disable all steppers but Z\n");

	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fclose(fp);
}

void  main(int argc, char** argv) {
	clock_t start, end;
	start = clock();//开始时间

	readfile(file_1);
	BoundingBox();//生成包围盒

	//Case1();
	Case2();
	Paths_Optimization();
	//Support();
	Contour();
	//GcodePrint();
	ContourGcodePrint();
	end = clock();//结束时间
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  //输出时间（单位：ｓ）

	//包围盒的路径
	//FILE* fp;
	//errno_t err;     //判断此文件流是否存在 存在返回1
	//err = fopen_s(&fp, "Box.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	//for (int i = 0; i < P.size(); i++) {
	//	for (int j = 0; j < P[i].size(); j++) {
	//		for (int m = 0; m < P[i][j].size(); m++) {
	//			if (m == 0)
	//				fprintf(fp, "%f %f %f %d\n", P[i][j][m].x, P[i][j][m].y, P[i][j][m].z, 0);
	//			else
	//				fprintf(fp, "%f %f %f %d\n", P[i][j][m].x, P[i][j][m].y, P[i][j][m].z, 1);
	//		}
	//	}
	//}
	//fclose(fp);

	//FILE* fp1;
	//err = fopen_s(&fp1, "contour.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	//for (int i = 0; i < shell1.size(); i++) {
	//	for (int j = 0; j < shell1[i].size(); j++) {
	//		for (int m = 0; m < shell1[i][j].size(); m++) {
	//			if (m == 0)
	//				fprintf(fp1, "%f %f %f %d\n", shell1[i][j][m].x, shell1[i][j][m].y, 2.0*(i+1), 0);
	//			else
	//				fprintf(fp1, "%f %f %f %d\n", shell1[i][j][m].x, shell1[i][j][m].y, 2.0 * (i + 1), 1);
	//		}
	//	}
	//}
	//fclose(fp1);

	//FILE* fp2;
	//err = fopen_s(&fp2, "paths.txt", "a"); //若return 1 , 则将指向这个文件的文件流给
	//for (int i = 0; i < paths.size(); i++) {
	//	for (int j = 0; j < paths[i].size(); j++) {
	//		for (int m = 0; m < paths[i][j].size(); m++) {
	//			if (m == 0)
	//				fprintf(fp2, "%f %f %f %d\n", paths[i][j][m].x, paths[i][j][m].y, paths[i][j][m].z, 0);
	//			else
	//				fprintf(fp2, "%f %f %f %d\n", paths[i][j][m].x, paths[i][j][m].y, paths[i][j][m].z, 1);
	//		}
	//	}
	//}
	//fclose(fp2);

}
