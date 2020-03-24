#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>
#include <windows.h>
#include <glut.h>
#include "bitmap_image.hpp"



using namespace std;

#define pi (2*acos(0.0))
#define epsilon 0.001

int currLevelOfRecursion = 0, maxLevelOfRecursion = 10;

class Vect
{
public :

	double x,y,z;

	    // constructs a vector with given components
    Vect(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vect() {
        this->x = 1;
        this->y = 1;
        this->z = 1;
    }

    void scaleBy (double factor){
	    x *= factor; y *= factor; z *= factor;
	}

	Vect getScaledVect (double factor){
	    Vect newVect;
	    newVect.x = x*factor; newVect.y = y*factor; newVect.z = z*factor;
	    return newVect;
	}

	void addTo (Vect a, double scalingFactor){
	    x += a.x*scalingFactor; y += a.y*scalingFactor; z += a.z*scalingFactor;
	}

    double getModulus(){
        return sqrt(x*x + y*y + z*z);
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    // add two vectors
    Vect operator+(const Vect& v)
    {
        Vect v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Vect operator-(const Vect& v)
    {
        Vect v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    // scale a vector with a given coefficient
    Vect operator* (double m)
    {
        Vect v(x*m, y*m, z*m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vect a, Vect b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    // get the cross product of two vectors
    static Vect cross(Vect a, Vect b)
    {
        Vect v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    // print a vector. only for testing purposes.
    void print ()
    {
        cout << "Vector : " <<  x << " " << y << " " << z << endl;
    }
};

int INVALID_POINT = 0, VALID_POINT = 1;

class Point
{
public:
	double x,y,z;
	int type;

    Point(double _x, double _y, double _z){
        x = _x; y = _y; z = _z;
        type = VALID_POINT;
    }

    Point(){
        x=1; y=1; z=1;
        type = VALID_POINT;
    }

    Point(Vect a){
        x=a.x; y=a.y; z=a.z;
    }

    void addTo (Vect a, double scalingFactor){
	    x += a.x*scalingFactor; y += a.y*scalingFactor; z += a.z*scalingFactor;
	}

    Point operator+(const Vect& v)
    {
        Point v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Point operator-(const Vect& v)
    {
        Point v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

/*    Vect operator+(const Point& v)
    {
        Vect v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }
*/
    // subtract one vector from another
    Vect operator-(const Point& v)
    {
        Vect v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    double getDistanceFrom(Point p){
        Vect dis = *this - p;
        return dis.getModulus();
    }

    void print ()
    {
        cout << "Point : " <<  x << " " << y << " " << z << endl;
    }
};



class Color {
public:
    double r, g, b;
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color() {
        r=0; g=0; b=0;
    }

    Color operator* (double m)
    {
        Color v(r*m, g*m, b*m);
        return v;
    }

    Color operator+(const Color& v)
    {
        Color v1(r+v.r, g+v.g, b+v.b);
        if (v1.r > 1) v1.r = 1;
        if (v1.g > 1) v1.g = 1;
        if (v1.b > 1) v1.b = 1;
        return v1;
    }

    void print(){
        cout << "Color : " << r << " " << g << " " << b << endl;
    }
};


int NORMAL_LIGHT = 0, SPOT_LIGHT = 1;

class Light{
public:
    Point position;
    int type;
    double cutOffAngle = 360;
    Vect direction;
    double falloff = .00001;


    Light(){
        position = Point(100,-200,100);
        type = NORMAL_LIGHT;
    }

    Light (Point pos, int ty, double fall){
        position = pos;
        type = ty;
        falloff = fall;
    }

    void setDirection(Vect v){
        direction = v;
        direction.normalize();
    }

    double getSourceIntensityAtPoint(Point p){
        double dist = position.getDistanceFrom(p);
        double i = exp(-falloff*dist*dist);
        if(type==SPOT_LIGHT){
            Vect objDir = (p - position); objDir.normalize();
            double angle = acos(Vect::dot(objDir, direction))*180/pi;
            if (abs(angle)>cutOffAngle) i = 0;
        }
        return i;
    }
};


Point cameraPos;
Vect l,r,u;
double fovX=80, fovY=80, aspectRatio=1, nearDistance=1, farDistance=1000;
double width, height;
int screenResolution;
int screenH = 768, screenW = 1366;

double cameraSpeed = 0.8;
double A = pi/60; double cosA = cos(A); double sinA = sin(A);

int slices = 25;

Color castRay (Point r0, Vect rd);


void initializeCamera (){
    cameraPos.x = -40; cameraPos.y = -170; cameraPos.z = 20;
    l.x = 0; l.y = 1; l.z = 0;
    r.x = 1; r.y = 0; r.z = 0;
    u.x = 0; u.y = 0; u.z = 1;

    height = 2 * nearDistance * tan(fovY * pi / 360.0);
    width = 2 * nearDistance * tan(fovY * aspectRatio * pi / 360.0);
}

void printLRU (){
    printf("\nl = (%f, %f, %f)\n", l.x,l.y,l.z);
    printf("r = (%f, %f, %f)\n", r.x,r.y,r.z);
    printf("u = (%f, %f, %f)\n", u.x,u.y,u.z);
}


// draw functions

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f(0, 0, 0);
		glVertex3f(0, a, 0);
		glVertex3f(a, a, 0);
		glVertex3f(a, 0, 0);
	}glEnd();
}


void drawSphere(double radius, Point center)
{
    Point c (center.x, center.y, center.z);
	struct Point points[100][100];
	int i,j;
	double h,r;
	int slices = 40, stacks = 40;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=-h;
		}
	}
	glPushMatrix();
	{

        glTranslatef(c.x, c.y, c.z);
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }
	}
	glPopMatrix();
}



class Object;
class Light;

vector<Object *> sceneObjects;
vector<Light *> sceneLights;


class Object{
public:
    double a, s, r, d, k;
    Color objColor;
    string type;

    Object(){}

    virtual double getIntersectingPoint(Point r0, Vect rd){}
    virtual Color getRawColorAtPoint(Point p){}
    virtual Vect getNormalAtGivenPoint (Point p){}
    virtual void drawObject(){}

    void setMaterialProperties(double _a, double _d, double _s, double _r, double _k){
        a = _a; d = _d; s = _s; r = _r; k = _k;
    }

    // ip = intersectionPoint
    Color getColorIntensity(Point r0, Vect rd, Point ip, Color currColor){
        double intensity = a;
        for (int l=0; l<sceneLights.size(); l++){
            Vect toLight = (sceneLights[l]->position - ip); toLight.normalize();
            Point startPoint = ip + toLight*epsilon;

            bool isIlluminated = true;
            for (int i=0; i<sceneObjects.size(); i++){
                double t = sceneObjects[i]->getIntersectingPoint(startPoint, toLight);
                if(t>0) {
                    isIlluminated = false; break;
                }
            }

            if (isIlluminated){
                double sourceScalingFactor = sceneLights[l]->getSourceIntensityAtPoint(ip);

                Vect N = getNormalAtGivenPoint(ip);
                double lambert = d * sourceScalingFactor * Vect::dot(toLight, N);

                Vect R = rd - N * (Vect::dot(rd, N) * 2); R.normalize();
                double phong = s * sourceScalingFactor * pow(Vect::dot(R, toLight), k);

                if (lambert < 0) lambert = 0; if (lambert > 1) lambert = 1;
                if (phong < 0) phong = 0; if (phong > 1) phong = 1;

                intensity += (phong + lambert);
                if (intensity>1) intensity = 1;
            }
        }
        currColor = currColor * intensity;

        if (currLevelOfRecursion == maxLevelOfRecursion) {
                currLevelOfRecursion = 0;
                return currColor;
        }
        else {
            currLevelOfRecursion++;
            Vect N = getNormalAtGivenPoint(ip);
            Vect reflect = rd - N * (Vect::dot(rd, N) * 2); reflect.normalize();
            Point newStart = ip + reflect*epsilon;
            return currColor + (castRay(newStart, reflect) * r);
        }
    }
};


class CheckerBoard : public Object{
public:
    int blkSize; double inf;
    CheckerBoard(int _blkSize, double _inf){
        blkSize = _blkSize; inf = _inf;
        a = 0.2; d = 0.2; s = 0.3; r = 0.3; k = 30;
        type = "Checkerboard";
    }


    void drawObject(){
        int blockCount = 0;
        for (int i=-inf; i<=inf; i+=blkSize){
            for (int j=-inf; j<=inf; j+=blkSize){
                if(blockCount%2) glColor3f(1,1,1);
                else glColor3f(0,0,0);
               // if ((i+inf)/blkSize == 20 && (j+inf)/blkSize == 20)
                 //   glColor3f(1,0,0);
                glPushMatrix();
                {
                    glTranslatef(i, j, 0);
                    drawSquare(blkSize);
                }
                glPopMatrix();
                blockCount+=1;
            }
        }
    }


    double getIntersectingPoint(Point r0, Vect rd){
        if (abs(rd.z)  < .00001 ){
            return -1;
        }
        return - r0.z / rd.z;
    }

    Color getRawColorAtPoint(Point p){
        int i = (p.x - (-inf))/blkSize;
        int j = (p.y - (-inf))/blkSize;

       // if(i==20 && j==20) return Color(1,0,0);
        if ((i+j)%2 == 0) return Color();
        else return Color(1,1,1);
    }

    Vect getNormalAtGivenPoint (Point p){
        return Vect(0,0,1);
    }
};


class Sphere : public Object{
public:
    double radius;
    Point c;
    Sphere(Point _c, double r){
        radius = r; c = _c;
        type = "Sphere";
    }


    void drawObject(){
        glColor3f(objColor.r, objColor.g, objColor.b);
        drawSphere(radius, c);
    }


    double getIntersectingPoint(Point r0, Vect rd){
        Vect R0 = r0 - c;
        double a = 1;
        double b = 2 * Vect::dot(rd, R0);
        double C = Vect::dot(R0, R0) - radius*radius;

        double d = b*b - 4*a*C;
        if (d<0) return -1;

        double t1 = (-b + sqrt(d)) / 2.0;
        double t2 = (-b - sqrt(d)) / 2.0;
        if (t1 > t2) return t2; else return t1;
    }

    Color getRawColorAtPoint(Point p){
        return objColor;
    }

    Vect getNormalAtGivenPoint (Point p){
        Vect N = (p-c); N.normalize();
        return N;
    }
};


class Triangle : public Object{
public:
    Point a, b, c;
    Triangle(Point _a, Point _b, Point _c){
        a = _a; b = _b; c = _c;
        type = "Triangle";
    }

    void drawObject(){
        glColor3f(objColor.r, objColor.g, objColor.b);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    static double getDetermiant (Vect a1, Vect a2, Vect a3){
        double det = 0;
        det += a1.x * (a2.y*a3.z - a3.y*a2.z);
        det -= a2.x * (a1.y*a3.z - a3.y*a1.z);
        det += a3.x * (a1.y*a2.z - a2.y*a1.z);
        return det;
    }

    double getIntersectingPoint(Point r0, Vect rd){
        double detA = Triangle::getDetermiant(a-b, a-c, rd);
        if(abs(detA) < 0.00001) return -1;

        double beta = Triangle::getDetermiant(a-r0, a-c, rd)/detA;
        double gamma = Triangle::getDetermiant(a-b, a-r0, rd)/detA;
        double t = Triangle::getDetermiant(a-b, a-c, a-r0)/detA;

        if (beta + gamma < 1 && beta>0 && gamma>0 && t>=0) return t;
        return -1;
    }

    Color getRawColorAtPoint(Point p){
        return objColor;
    }

    Vect getNormalAtGivenPoint (Point p){
        Vect pa = a - p;
        Vect pb = b - p;
        Vect N = Vect::cross(pa, pb); N.normalize();
        return N;
    }
};



Color castRay (Point r0, Vect rd){
    double tMin = 99999; int objId = -1;
    for (int k = 0; k<sceneObjects.size(); k++){
            double t = sceneObjects[k]->getIntersectingPoint(r0, rd);
            if (t<0) continue;
            else if (cameraPos.getDistanceFrom(r0 + rd*t) > farDistance) continue;
            if (t < tMin){
                objId = k;
                tMin = t;
            }
    }
    Color pixelColor = Color();
    if (objId != -1){
        pixelColor = sceneObjects[objId]->getColorIntensity(r0, rd, (r0 + rd*tMin), sceneObjects[objId]->getRawColorAtPoint((r0 + rd*tMin)));
    } else {
        currLevelOfRecursion = 0;
    }
    return pixelColor;
}



void saveScreenshot (){
    cout << "Capturing..." << endl;
    Point sceneMidPoint = cameraPos + l*nearDistance;
    Point lowLeftCorner = sceneMidPoint - u*(height/2.0) - r*(width/2.0);

    double pw = width/screenW, ph = height/screenH;

    Color pixelColor;
    int total = screenW * screenH;
    int pc = 0;
    bitmap_image image(screenW, screenH);

    for (int i=0; i<screenW; i++){
        for (int j=0; j<screenH; j++){
            Point r0 (lowLeftCorner + r*(i*pw) + u*(j*ph));
            Vect rd = r0 - cameraPos;
            rd.normalize();
            pixelColor = castRay(r0, rd);
            image.set_pixel(i, screenH-j-1, pixelColor.r*255, pixelColor.g*255, pixelColor.b*255);
            pc++;
            if (pc% (total/5) == 0) {
                int percentage = int (pc*100.0/total); cout << percentage << "% done..." << endl;
            }
        }
    }
    image.save_image("out.bmp");
    cout << "Finished. Saved to out.bmp." << endl;
}


//listeners

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case '1':{
            Vect lCosA = l.getScaledVect(cosA);
            Vect lSinA = l.getScaledVect(sinA);
            Vect rCosA = r.getScaledVect(cosA);
            Vect rSinA = r.getScaledVect(sinA);

            l.x = lCosA.x + rSinA.x;
            l.y = lCosA.y + rSinA.y;
            l.z = lCosA.z + rSinA.z;

            r.x = rCosA.x - lSinA.x;
            r.y = rCosA.y - lSinA.y;
            r.z = rCosA.z - lSinA.z;
            printLRU();
			break;
		}
		case '2':{
            Vect lCosA = l.getScaledVect(cosA);
            Vect lSinA = l.getScaledVect(sinA);
            Vect rCosA = r.getScaledVect(cosA);
            Vect rSinA = r.getScaledVect(sinA);

            l.x = lCosA.x - rSinA.x;
            l.y = lCosA.y - rSinA.y;
            l.z = lCosA.z - rSinA.z;

            r.x = rCosA.x + lSinA.x;
            r.y = rCosA.y + lSinA.y;
            r.z = rCosA.z + lSinA.z;
            printLRU();
			break;
		}
        case '3':{
            Vect lCosA = l.getScaledVect(cosA);
            Vect lSinA = l.getScaledVect(sinA);
            Vect uCosA = u.getScaledVect(cosA);
            Vect uSinA = u.getScaledVect(sinA);

            l.x = lCosA.x + uSinA.x;
            l.y = lCosA.y + uSinA.y;
            l.z = lCosA.z + uSinA.z;

            u.x = uCosA.x - lSinA.x;
            u.y = uCosA.y - lSinA.y;
            u.z = uCosA.z - lSinA.z;
            printLRU();
			break;
		}
		case '4':{
            Vect lCosA = l.getScaledVect(cosA);
            Vect lSinA = l.getScaledVect(sinA);
            Vect uCosA = u.getScaledVect(cosA);
            Vect uSinA = u.getScaledVect(sinA);

            l.x = lCosA.x - uSinA.x;
            l.y = lCosA.y - uSinA.y;
            l.z = lCosA.z - uSinA.z;

            u.x = uCosA.x + lSinA.x;
            u.y = uCosA.y + lSinA.y;
            u.z = uCosA.z + lSinA.z;
            printLRU();
			break;
		}
		case '6':{
            Vect rCosA = r.getScaledVect(cosA);
            Vect rSinA = r.getScaledVect(sinA);
            Vect uCosA = u.getScaledVect(cosA);
            Vect uSinA = u.getScaledVect(sinA);

            r.x = rCosA.x + uSinA.x;
            r.y = rCosA.y + uSinA.y;
            r.z = rCosA.z + uSinA.z;

            u.x = uCosA.x - rSinA.x;
            u.y = uCosA.y - rSinA.y;
            u.z = uCosA.z - rSinA.z;
            printLRU();
			break;
		}
		case '5':{
            Vect rCosA = r.getScaledVect(cosA);
            Vect rSinA = r.getScaledVect(sinA);
            Vect uCosA = u.getScaledVect(cosA);
            Vect uSinA = u.getScaledVect(sinA);

            r.x = rCosA.x - uSinA.x;
            r.y = rCosA.y - uSinA.y;
            r.z = rCosA.z - uSinA.z;

            u.x = uCosA.x + rSinA.x;
            u.y = uCosA.y + rSinA.y;
            u.z = uCosA.z + rSinA.z;
            printLRU();
			break;
		}
		case '0':{
            saveScreenshot();
		}
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraPos.addTo(l, -cameraSpeed);
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraPos.addTo(l, cameraSpeed);
			break;

		case GLUT_KEY_RIGHT:
			cameraPos.addTo(r, cameraSpeed);
			break;
		case GLUT_KEY_LEFT:
			cameraPos.addTo(r, -cameraSpeed);
			break;

		case GLUT_KEY_PAGE_UP:
		    cameraPos.addTo(u, cameraSpeed);
			break;
		case GLUT_KEY_PAGE_DOWN:
		    cameraPos.addTo(u, -cameraSpeed);
			break;
		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
			}
			break;
		default:
			break;
	}
}


void display(){
	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);
	//initialize the matrix
	glLoadIdentity();

    gluLookAt(cameraPos.x, cameraPos.y, cameraPos.z,
              cameraPos.x + l.x, cameraPos.y + l.y, cameraPos.z + l.z,
              u.x, u.y, u.z);
	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);
	/****************************
	/ Add your objects from here
	****************************/
	//add objects
    for (int i=0; i<sceneObjects.size(); i++){
        sceneObjects[i]->drawObject();
    }
    glColor3f(0.8,0.8,0.8);
    for (int i=0; i<sceneLights.size(); i++){
        drawSphere(5, sceneLights[i]->position);
    }
    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//clear the screen
	glClearColor(0,0,0,0);
	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fovY, aspectRatio, nearDistance,	farDistance);
    initializeCamera();
}

void setTheStageFromFileDescription(){
    ifstream descript;
    descript.open("description.txt");
    descript >> nearDistance >> farDistance >> fovY >> aspectRatio;
    descript >> maxLevelOfRecursion;
    descript >> screenResolution;
    aspectRatio = 1366.0/768.0;

    double a, d, s, r, k, x, y, z, R, G, B, fallOffParam, pyrH, pyrW;
    //checker board
    int blockSize = 0;
    descript >> blockSize;
    sceneObjects.push_back(new CheckerBoard(blockSize, blockSize*20));
    descript >> sceneObjects[0]->a >> sceneObjects[0]->d >> sceneObjects[0]->r;
    sceneObjects[0]->s = 0.2;

    int objCount = 0;
    descript >> objCount;
    string type;
    for(int i=1; i<objCount+1; i++)
    {
        descript >> type;

        if(type == "sphere")
        {
            descript >> x >> y >> z;
            descript >> r;
            sceneObjects.push_back(new Sphere(Point(x, y, z), r));

            descript >> R >> G >> B;
            sceneObjects[i]->objColor = Color(R, G, B);

            descript >> a >> d >> s >> r;
            descript >> k;
            sceneObjects[i]->setMaterialProperties(a, d, s, r, k);
        }

        else if(type == "pyramid")
        {
            descript >> x >> y >> z;
            descript >> pyrW >> pyrH;
            descript >> R >> G >> B;
            descript >> a >> d >> s >> r;
            descript >> k;

            Point p1(x, y, z), p2(x-pyrW, y, z), p3(x-pyrW, y+pyrW, z), p4(x, y+pyrW, z), p5(x-pyrW/2.0, y+pyrW/2.0, z+pyrH);

            sceneObjects.push_back(new Triangle(p1, p2, p3));
            sceneObjects.push_back(new Triangle(p1, p3, p4));
            sceneObjects.push_back(new Triangle(p1, p2, p5));
            sceneObjects.push_back(new Triangle(p2, p3, p5));
            sceneObjects.push_back(new Triangle(p3, p4, p5));
            sceneObjects.push_back(new Triangle(p4, p1, p5));

            for (int j = 0; j<6; j++){
                sceneObjects[i]->setMaterialProperties(a, d, s, r, k);
                sceneObjects[i]->objColor = Color(R, G, B);
                i++; objCount++;
            }
            i--; objCount--;
        }
    }
    cout << "Objects loaded into scene." << endl;

    int lightCount = 0;
    descript >> lightCount;
    for(int i=0; i<lightCount; i++)
    {
        descript >> x >> y >> z >> fallOffParam;
        sceneLights.push_back(new Light(Point(x, y, z), NORMAL_LIGHT, fallOffParam));
    }

    descript >> lightCount;
    for(int j=0; j<lightCount; j++)
    {
        descript >> x >> y >> z >> fallOffParam;
        Light* spot = new Light(Point(x, y, z), SPOT_LIGHT, fallOffParam);

        double cAngle;
        descript >> x >> y >> z >> cAngle;

        Vect direction = Point(x, y, z) - spot->position;
        spot->setDirection(direction);
        spot->cutOffAngle = cAngle;
        sceneLights.push_back(spot);
    }
    cout << "Lights loaded into scene." << endl;
    descript.close();
    cout << "Finished reading input. File closed." << endl;
}


int main(int argc, char **argv){
    setTheStageFromFileDescription();

	glutInit(&argc,argv);
	glutInitWindowSize(screenW/2, screenH/2);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); //Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracing");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

    sceneLights.clear();
    sceneObjects.clear();

	return 0;
}
