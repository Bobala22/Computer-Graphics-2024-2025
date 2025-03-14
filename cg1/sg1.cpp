/*
  This program plots different 2D functions.
*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>
// #include "glut.h" //MSVC local library install
#include <GL/glut.h> //system-wide install (or compiler default path)
#include <vector>

double circle = atan(1) * 8;
double halfCircle = atan(1) * 4;
double tau = circle;    // 2 * PI = TAU
double pi = halfCircle; // TAU / 2 = PI

// How often should the drawing algorithm sample the function.
double step = 0.05;

int defaultW = 1000, defaultH = 1000;

unsigned char prevKey;

/*
   Nicomedes' Conchoid
   $x = a + b \cdot cos(t), y = a \cdot tg(t) + b \cdot sin(t)$. or
   $x = a - b \cdot cos(t), y = a \cdot tg(t) - b \cdot sin(t)$. where
   $t \in (-\pi / 2, \pi / 2)$
*/
void Display1()
{
  double xmax, ymax, xmin, ymin;
  /*
    Nicomedes' Conchoid is a family of functions. These two parameters,
    a and b, choose a specific 2D function from that family.
    It has two branches, so, for each y, we have two x values.
    Therefore, we're actually drawing two lines, not one.
   */
  double a = 1, b = 2;

  /*
    We shadow the global variable with a different step size,
    as values different from 0.05
    would yield a significantly different plot.
  */
  double step = 0.05;

  /*
    First, we compute the points of the function, so we can determine
    the maximal extend of the drawing.
   */
  xmax = a - b - 1;
  xmin = a + b + 1;
  ymax = ymin = 0;
  for (double t = -pi / 2 + step; t < pi / 2; t += step)
  {
    double x1, y1, x2, y2;
    x1 = a + b * cos(t);
    xmax = (xmax < x1) ? x1 : xmax;
    xmin = (xmin > x1) ? x1 : xmin;

    x2 = a - b * cos(t);
    xmax = (xmax < x2) ? x2 : xmax;
    xmin = (xmin > x2) ? x2 : xmin;

    y1 = a * tan(t) + b * sin(t);
    ymax = (ymax < y1) ? y1 : ymax;
    ymin = (ymin > y1) ? y1 : ymin;

    y2 = a * tan(t) - b * sin(t);
    ymax = (ymax < y2) ? y2 : ymax;
    ymin = (ymin > y2) ? y2 : ymin;
  }
  // We care about the maximal extent on each axis (from the origin).
  xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
  ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

  /*
    Since we have, for x and for y, the maximal absolute values,
    dividing the coordinates of the points by these values will
    ensure we draw the whole function inside [-1, 1]^2, the default
    OpenGL screen.
   */

  glColor3f(1, 0.1, 0.1);
  glBegin(GL_LINE_STRIP);
  for (double t = -pi / 2 + step; t < pi / 2; t += step)
  {
    double x1, y1, x2, y2; // You might get some warnings in your IDE. Why?
    x1 = (a + b * cos(t)) / xmax;
    x2 = (a - b * cos(t)) / xmax;
    y1 = (a * tan(t) + b * sin(t)) / ymax;
    y2 = (a * tan(t) - b * sin(t)) / ymax;

    glVertex2d(x1, y1);
  }
  glEnd();

  glBegin(GL_LINE_STRIP);
  for (double t = -pi / 2 + step; t < pi / 2; t += step)
  {
    double x1, y1, x2, y2;
    x1 = (a + b * cos(t)) / xmax;
    x2 = (a - b * cos(t)) / xmax;
    y1 = (a * tan(t) + b * sin(t)) / ymax;
    y2 = (a * tan(t) - b * sin(t)) / ymax;

    glVertex2d(x2, y2);
  }
  glEnd();
}

// $f(x) = \left| sin(x) \right| \cdot e^{-sin(x)}, x \in \left[ 0, 8 \cdot \pi \right]$,
void Display2()
{
  /*
    We can determine how far the function extends
    (and thus the needed scaling factors)
    by looking at the function and doing a bit of Calculus.
  */
  double xmax = 8 * pi;
  double ymax = exp(1.1); // Why 1.1?

  glColor3f(1, 0.1, 0.1);
  glBegin(GL_LINE_STRIP);
  for (double x = 0; x < xmax; x += step)
  {
    double x1, y1;
    x1 = x / xmax;
    y1 = (fabs(sin(x)) * exp(-sin(x))) / ymax;
    glVertex2d(x1, y1);
  }
  glEnd();
}

/* 1)
   \( f(x) =
     \left\{
       \begin{array}{cl}
         1              & x = 0   \\
         \frac{d(x)}{x} & x \gt 0 \\
       \end{array}
     \right.
   \)
 */
void Display3()
{

  auto d = [](double x)
  {
    double nearest_int = round(x);
    return fabs(x - nearest_int);
  };

  double xmax = 20.0;
  double ymax = 1.0;

  glColor3f(1, 0.1, 0.1);
  glBegin(GL_LINE_STRIP);

  glVertex2d(0, 0.9);
  glVertex2d(0.02, 0.9);

  for (double x = 1; x <= 100; x += step)
  {
    double x1, y1;
    x1 = x / xmax;
    y1 = (d(x) / x) / ymax;
    glVertex2d(x1, y1);
  }

  glEnd();
}

// 3) function arguments e.g.: f(a, b, t), where a and b are function family parameters, and the is the driving variables.
void plot(double (*x)(double, double, double), double (*y)(double, double, double), double a, double b, double intervalStart, double intervalEnd, double step = 0.01, double scaleX = 1, double scaleY = 1, GLint primitive = GL_LINE_STRIP)
{
  glColor3f(1, 0.1, 0.1);
  glBegin(GL_LINE_STRIP);

  for (double t = intervalStart; t <= intervalEnd; t += step)
  {
    double x1 = x(a, b, t) / scaleX;
    double y1 = y(a, b, t) / scaleY;
    glVertex2d(x1, y1);
  }

  glEnd();
  glFlush();
}

// Function definitions for each curve

// **1) Circle Conchoid (Lima�on, Pascal's Snail)**
double x_circle_conchoid(double a, double b, double t)
{
  return 2 * (a * cos(t) + b) * cos(t);
}

double y_circle_conchoid(double a, double b, double t)
{
  return 2 * (a * cos(t) + b) * sin(t);
}

// **2) Cycloid**
double x_cycloid(double a, double b, double t)
{
  return a * t - b * sin(t);
}

double y_cycloid(double a, double b, double t)
{
  return a - b * cos(t);
}

// **3) Epicycloid**
double x_epicycloid(double a, double b, double t)
{
  return (a + b) * cos((b / a) * t) - b * cos(t + (b / a) * t);
}

double y_epicycloid(double a, double b, double t)
{
  return (a + b) * sin((b / a) * t) - b * sin(t + (b / a) * t);
}

// **4) Hypocycloid**
double x_hypocycloid(double a, double b, double t)
{
  return (a - b) * cos((b / a) * t) - b * cos(t - (b / a) * t);
}

double y_hypocycloid(double a, double b, double t)
{
  return (a - b) * sin((b / a) * t) - b * sin(t - (b / a) * t);
}

// 3
/*
  2) Circle Concoid (Lima�on, Pascal's Snail):
  \(x = 2 \cdot (a \cdot cos(t) + b) \cdot cos(t), \; y = 2 \cdot (a \cdot cos(t) + b) \cdot sin(t), \; t \in (-\pi, \pi)\) .
  For this plot, \(a = 0.3, \; b = 0.2\) .
*/
void Display4()
{
  plot(x_circle_conchoid, y_circle_conchoid, 0.3, 0.2, -pi, pi, 0.01, 1.5, 0.7);
}

/*
  2) Cicloid:
  \(x = a \cdot t - b \cdot sin(t), \; y = a - b \cdot cos(t), \; t \in \mathbb{R} \) .
  For this plot, \(a = 0.1, \; b = 0.2\) .
*/
void Display5()
{
  plot(x_cycloid, y_cycloid, 0.1, 0.2, -3 * pi + 0.1, 3 * pi - 0.1, 0.01, 1.0, 1.2);
}

/*
  2) Epicicloid:
  \(x = (a + b) \cdot cos\left( \frac{b}{a} \cdot t \right) - b \cdot cos\left(t + \frac{b}{a}\cdot t \right) \)
  \(y = (a + b) \cdot sin\left( \frac{b}{a} \cdot t \right) - b \cdot sin\left(t + \frac{b}{a}\cdot t \right) \)
  \( t \in \left[ 0, 2\pi \right] \) .
  For this plot, \(a = 0.1, \; b = 0.3\) .
*/
void Display6()
{
  plot(x_epicycloid, y_epicycloid, 0.1, 0.3, 0, 2 * pi, 0.01, 0.8, 0.8);
}

/*
  2) Hipocicloid:
  \(x = (a - b) \cdot cos\left( \frac{b}{a} \cdot t \right) - b \cdot cos\left(t - \frac{b}{a}\cdot t \right) \)
  \(y = (a - b) \cdot sin\left( \frac{b}{a} \cdot t \right) - b \cdot sin\left(t - \frac{b}{a}\cdot t \right) \)
  \( t \in \left[ 0, 2\pi \right] \) .
  For this plot, \(a = 0.1, \; b = 0.3\) .
 */
void Display7()
{
  plot(x_hypocycloid, y_hypocycloid, 0.1, 0.3, 0, 2 * pi, 0.01, 0.6, 0.6);
}

/*
 4) Logarithmic spiral (in polar coordinates):
 \( r = a \cdot e^{1+t}, \; t \in (0, \infty) \) .
 For this plot, \(a = 0.02\) .
*/

double x_spiral(double a, double b, double t)
{
  double r = a * exp(1 + t);
  return r * cos(t);
}

double y_spiral(double a, double b, double t)
{
  double r = a * exp(1 + t);
  return r * sin(t);
}

void Display8()
{

  plot(x_spiral, y_spiral, 0.02, 0, -2, 5);
}

/*
  4) Sine polar plot flower:
  \( r = sin(a \cdot t), \; t \in (0, \infty)  \) .
  For this plot, \(a = 10\), and the number 'petals' is \( 2 \cdot a \). Think about why.
*/

double x_flower(double a, double b, double t)
{
  double r = sin(a * t);
  return r * cos(t);
}

double y_flower(double a, double b, double t)
{
  double r = sin(a * t);
  return r * sin(t);
}

void Display9()
{
  plot(x_flower, y_flower, 10, 0, 0, 9, 0.01, 1.1, 1.1);
}

/*
5) Longchamps' Trisectrix:
\(
x = \frac{a}{4 \cdot cos^2(t) - 3}, \;
y = \frac{a \cdot tg(t)}{4 \cdot cos^2(t) - 3}, \;
t \in (-\pi/2, \pi/2) \setminus \{ -\pi/6, \pi/6 \} \) .
For this plot, \(a = 0.2\) .
 */
void plotTrisectrix(double (*xFunc)(double, double, double),
                    double (*yFunc)(double, double, double),
                    double a,
                    double b,
                    double intervalStart = -M_PI / 2,
                    double intervalEnd = M_PI / 2,
                    double step = 0.03,
                    double scaleX = 1,
                    double scaleY = 1,
                    GLint /*primitive*/ = GL_LINE_STRIP) // primitive not used here
{
  // Local structure to hold vertices.
  struct Point2D
  {
    double x, y;
  };

  // Boundaries for clamping.
  const double xMin = -0.95;
  const double yMax = 0.95;

  // Offset to move the plot slightly down and to the RIGHT.
  const double offsetX = 0.05;
  const double offsetY = -0.05;

  // Container to collect curve vertices.
  std::vector<Point2D> curveVertices;

  // Generate vertices along the trisectrix.
  for (double t = intervalStart; t <= intervalEnd; t += step)
  {
    // Skip the excluded t values.
    if (fabs(t + M_PI / 6) < 1e-9 || fabs(t - M_PI / 6) < 1e-9)
      continue;

    // Compute raw values.
    double rawX = xFunc(a, b, t) * scaleX;
    double rawY = yFunc(a, b, t) * scaleY;

    // Clamp rawY to the top boundary.
    if (rawY > yMax)
      rawY = yMax;

    // Force x to xMin if rawX is too far left.
    double newX = (rawX < xMin) ? xMin : rawX;
    double newY = rawY;

    // Only consider points in the left (x <= 0) and top (y >= 0) quadrant.
    if (!(newX <= 0.0 && newY >= 0.0))
      continue;

    // Store vertex with translation.
    curveVertices.push_back({newX + offsetX, newY + offsetY});
  }

  // Need at least two vertices to form triangles.
  if (curveVertices.size() < 2)
    return;

  // Fixed corner at the intersection of left and top boundaries.
  Point2D corner = {xMin + offsetX, yMax + offsetY};

  // Draw triangles with alternating colors.
  int numTriangles = curveVertices.size() - 1;
  int mid = numTriangles / 2;

  glBegin(GL_TRIANGLES);
  for (int i = 0; i < numTriangles; i++)
  {
    // First half: alternate black and white;
    // Second half: alternate red and white.
    if (i < mid)
    {
      if (i % 2 == 0)
        glColor3f(0.0f, 0.0f, 0.0f); // black
      else
        glColor3f(1.0f, 1.0f, 1.0f); // white
    }
    else
    {
      if (i % 2 == 0)
        glColor3f(1.0f, 0.0f, 0.0f); // red
      else
        glColor3f(1.0f, 1.0f, 1.0f); // white
    }
    glVertex2d(corner.x, corner.y);
    glVertex2d(curveVertices[i].x, curveVertices[i].y);
    glVertex2d(curveVertices[i + 1].x, curveVertices[i + 1].y);
  }
  glEnd();

  // Draw the border curved line along the collected vertices.
  size_t midIndex = curveVertices.size() / 2;

  glLineWidth(2);

  // Draw first half in black.
  glColor3f(0.0f, 0.0f, 0.0f); // black
  glBegin(GL_LINE_STRIP);
  for (size_t i = 0; i < midIndex; i++)
  {
    glVertex2d(curveVertices[i].x, curveVertices[i].y);
  }
  glEnd();

  // Draw second half in red.
  glColor3f(1.0f, 0.0f, 0.0f); // red
  glBegin(GL_LINE_STRIP);
  for (size_t i = midIndex; i < curveVertices.size(); i++)
  {
    glVertex2d(curveVertices[i].x, curveVertices[i].y);
  }
  glEnd();

  glFlush();
}

double x_trisectrix(double a, double b, double t)
{
  return a / (4 * pow(cos(t), 2) - 3);
}

double y_trisectrix(double a, double b, double t)
{
  return (a * tan(t)) / (4 * pow(cos(t), 2) - 3);
}

void Display10()
{
  plotTrisectrix(x_trisectrix, y_trisectrix, 0.2, 0, -pi / 2, pi / 2, 0.03, 1.0, 1.0, GL_LINE_STRIP);
}

void init(void)
{
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glLineWidth(2);
  glPointSize(1);
  // glPolygonMode(GL_FRONT, GL_LINE);
  // Enabling blending and smoothing
  glEnable(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_NICEST, GL_POINT_SMOOTH_HINT);
  glHint(GL_NICEST, GL_LINE_SMOOTH_HINT);
  glHint(GL_NICEST, GL_POLYGON_SMOOTH_HINT);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
}

void Display(void)
{
  std::cout << ("Call Display") << std::endl;
  // Clear the buffer. See init();
  glClear(GL_COLOR_BUFFER_BIT);

  switch (prevKey)
  {
  case '1':
    Display1();
    break;
  case '2':
    Display2();
    break;
  case '3':
    Display3();
    break;
  case '4':
    Display4();
    break;
  case '5':
    Display5();
    break;
  case '6':
    Display6();
    break;
  case '7':
    Display7();
    break;
  case '8':
    Display8();
    break;
  case '9':
    Display9();
    break;
  case '0':
    Display10();
    break;
  default:
    break;
  }
  glFlush();
}

void Reshape(int w, int h)
{
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}
void KeyboardFunc(unsigned char key, int x, int y)
{
  prevKey = key;
  if (key == 27) // escape
    exit(0);
  // The proper way to ask glut to redraw the window.
  glutPostRedisplay();
}

/*
  Callback upon mouse press or release.
  The button can be:
  GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, GLUT_RIGHT_BUTTON
  (and further for mousewheel and other mouse buttons)
  The state can be either GLUT_DOWN or  GLUT_UP, for
  a pressed or released button.
  (x, y) are the coordinates of the mouse.
*/
void MouseFunc(int button, int state, int x, int y)
{
  std::cout << "Mouse button ";
  std::cout << ((button == GLUT_LEFT_BUTTON) ? "left" : ((button == GLUT_RIGHT_BUTTON) ? "right" : "middle")) << " ";
  std::cout << ((state == GLUT_DOWN) ? "pressed" : "released");
  std::cout << " at coordinates: " << x << " x " << y << std::endl;
}

int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitWindowSize(defaultW, defaultH);
  glutInitWindowPosition(-1, -1);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
  glutCreateWindow(argv[0]);
  init();
  glutReshapeFunc(Reshape);
  glutKeyboardFunc(KeyboardFunc);
  glutMouseFunc(MouseFunc);
  glutDisplayFunc(Display);
  // glutIdleFunc(Display);
  glutMainLoop();

  return 0;
}