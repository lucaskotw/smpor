#include "draw_layout.h"


/*******************
 * Basic Operation *
 *******************/
CoordType absolute(CoordType coordVal)
{
    if (coordVal < 0) return -coordVal;
    else return coordVal;
}

double get_drawing_width(std::vector< std::vector<CoordType> >& coord)
{
    double max_width = MAX_WIDTH_INIT_VAL;
    for (std::vector< std::vector<CoordType> >::iterator itr=coord.begin();\
        itr!=coord.end();\
        ++itr)
    {
        for (std::vector<CoordType>::iterator itc=(*itr).begin(); \
            itc!=(*itr).end(); \
            ++itc)
        {
            if ( absolute( *itc ) > max_width) max_width = *itc;
        }
    }
    return max_width;
}




/******************
 * Basic Property *
 ******************/
void set_points_attributes()
{
  glEnable(GL_POINT_SMOOTH);
  glPointSize(10);
}





/******************
 * Drawing Method *
 ******************/
void draw_vertices(std::vector< std::vector<CoordType> >& coord, \
    std::vector<PartType>& partition)
{
  double coord_max = get_drawing_width(coord);

  glBegin(GL_POINTS);

    for (int i=0; i<coord[0].size(); ++i) {
        switch(partition.at(i))
        {
            case 0:
                glColor3f(1.0f, 0.0f, 0.0f);  
                break;
            case 1:
                glColor3f(0.0f, 1.0f, 0.0f);  
                break;
            case 2:
                glColor3f(0.0f, 0.0f, 1.0f);  
                break;
            case 3:
                glColor3f(0.0f, 0.0f, 0.0f);  
                break;
            default:
                glColor3f(0.5f, 0.5f, 0.5f);  
                break;
        }

        glVertex2f(coord[i][0]/coord_max, coord[i][1]/coord_max);
    }

  glEnd();
}


/*
 * This is only for "straight" edge
 */
void draw_edges(std::vector< std::vector<VtxType> >& edges, \
    std::vector< std::vector<CoordType> >& coord)
{
  double coord_max = get_drawing_width(coord);
  glBegin(GL_LINES);

    glColor3f(0.5f, 0.5f, 0.5f);
    double pt1;
    double pt2;
    for (int i=0; i<edges[0].size(); ++i) {
      pt1 = edges[i][0];
      pt2 = edges[i][1];
      glVertex2f(coord[pt1][0]/coord_max, coord[pt1][1]/coord_max);
      glVertex2f(coord[pt2][0]/coord_max, coord[pt2][1]/coord_max);
      // std::cout << "edges size=" << i << std::endl;
    }
    

  glEnd();
}


void draw_layout(std::vector< std::vector<VtxType> >& edges, \
    std::vector< std::vector<CoordType> >& coord,
    std::vector<PartType>& partition)
{


  GLFWwindow* window;

  /* Initialize the library */
  if (!glfwInit())
      exit(1);
      // return -1;

  /* Create a windowed mode window and its OpenGL context */
  window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "2D Stress Layout", NULL, NULL);
  if (!window)
  {
      glfwTerminate();
      exit(1);
      // return -1;
  }

  /* Make the window's context current */
  glfwMakeContextCurrent(window);


  /* Loop until the user closes the window */
  while (!glfwWindowShouldClose(window))
  {
    /* Render here */

    // glfwGetFramebufferSize(window, &width, &height);
    // glViewport(0, 0, width, height);
    set_points_attributes();

    // set frame background color
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear(GL_COLOR_BUFFER_BIT);


    draw_vertices(coord, partition);
    draw_edges(edges, coord);

    


    /* Swap front and back buffers */
    glfwSwapBuffers(window);

    /* Poll for and process events */
    glfwPollEvents();
  }

  glfwTerminate();
}