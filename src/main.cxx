// generated by Fast Light User Interface Designer (fluid) version 1.0110

#include "main.h"

static void timeout_cb() {
  glWindow->rotateY++;
	if (glWindow->rotateY >= 360)
	{
		glWindow->rotateY -= 360;
	}
	
	glWindow->redraw();
Fl::repeat_timeout(TIMESTEP, (Fl_Timeout_Handler)timeout_cb);
}

Gl_Window *glWindow=(Gl_Window *)0;

static void cb_Select(Fl_Button*, void*) {
  glWindow->mesh->LoadObjFile(fl_file_chooser("Choose an Object File", "*.obj", NULL));
//disSlider->value(0);
glWindow->isNewWeaving = true;
glWindow->isNewBaseMesh = true;
glWindow->redraw();
} 

static void cb_Background(Fl_Light_Button* o, void*) {
  if (o->value())
	  glClearColor(1.0, 1.0, 1.0, 1.0);
  else
	  glClearColor(0.0, 0.0, 0.0, 0.0);
  
  glWindow->redraw();
}

static void cb_Base(Fl_Light_Button* o, void*) {
  glWindow->isBaseMesh = bool(o->value());
glWindow->redraw();
}

static void cb_Texture(Fl_Button*, void*) {
  glWindow->setTexture( fl_file_chooser("Choose an Texture File", NULL, NULL), 1 ); // absolute value
glWindow->redraw();
}

static void cb_Environment(Fl_Button*, void*) {
  glWindow->setTexture( fl_file_chooser("Choose an Texture File", NULL, NULL), 0 ); // absolute value
glWindow->redraw();
}

static void cb_6(Fl_Light_Button* o, void*) {
  if(o->value()) // down - begin to play
    {
        Fl::add_timeout(TIMESTEP, (Fl_Timeout_Handler)timeout_cb);
        o->label("   @+3||"); // change to ||
    }
    else // up - begin to pause
    {
        Fl::remove_timeout((Fl_Timeout_Handler)timeout_cb);
        o->label("   @+6>"); // change to |>
    };
}

static void cb_Displace(Fl_Value_Slider* o, void*) {
  glWindow->mesh->SetDiplaceFactor(o->value());
glWindow->isNewWeaving = true;
glWindow->redraw();
}

static void cb_Width(Fl_Value_Slider* o, void*) {
  glWindow->mesh->SetWidth(o->value());
glWindow->isNewWeaving = true;
glWindow->redraw();
}

static void cb_Curvature(Fl_Value_Slider* o, void*) {
  glWindow->mesh->SetCurvature(o->value());
glWindow->isNewWeaving = true;
glWindow->redraw();
}

static void cb_center(Fl_Value_Slider* o, void*) {
  glWindow->mesh->SetCenterWeight(o->value());
glWindow->isNewWeaving = true;
glWindow->redraw();
}

static void cb_flatten(Fl_Button*, void*) {
  glWindow->mesh->Flatten();
glWindow->isNewBaseMesh = true;
glWindow->redraw();
}

static void cb_write_file(Fl_Button*, void*) {
  glWindow->mesh->Placing2D();
glWindow->mesh->OutputVertexPaperStrips();
glWindow->isNewBaseMesh = true;
glWindow->redraw();
}

static void cb_hole(Fl_Value_Input* o, void*) {
  glWindow->mesh->holeRadius = o->value();
}

static void cb_paper(Fl_Value_Input* o, void*) {
  glWindow->mesh->paperWidth = o->value();
}

static void cb_paper1(Fl_Value_Input* o, void*) {
  glWindow->mesh->paperHeight = o->value();
}

static void cb_piece(Fl_Value_Input* o, void*) {
  glWindow->mesh->boundingBox_size = o->value();
}

int main(int argc, char **argv) {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = new Fl_Double_Window(1235, 799, "Weaving");
    w = o;
    { glWindow = new Gl_Window(0, 0, 1000, 800, "label");
      glWindow->box(FL_PLASTIC_DOWN_FRAME);
      glWindow->color((Fl_Color)FL_BACKGROUND_COLOR);
      glWindow->selection_color((Fl_Color)FL_BACKGROUND_COLOR);
      glWindow->labeltype(FL_NORMAL_LABEL);
      glWindow->labelfont(0);
      glWindow->labelsize(14);
      glWindow->labelcolor((Fl_Color)FL_FOREGROUND_COLOR);
      glWindow->align(FL_ALIGN_CENTER);
      glWindow->when(FL_WHEN_RELEASE);
    } // Gl_Window* glWindow
    { Fl_Button* o = new Fl_Button(1005, 15, 95, 40, "Select Object");
      o->box(FL_GTK_UP_BOX);
      o->down_box(FL_GTK_DOWN_BOX);
      o->callback((Fl_Callback*)cb_Select);
    } // Fl_Button* o
    { Fl_Light_Button* o = new Fl_Light_Button(1005, 80, 95, 40, "Background");
      o->box(FL_GTK_UP_BOX);
      o->down_box(FL_GTK_DOWN_BOX);
      o->value(1);
      o->selection_color((Fl_Color)1);
      o->callback((Fl_Callback*)cb_Background);
    } // Fl_Light_Button* o
    { Fl_Light_Button* o = new Fl_Light_Button(1005, 130, 95, 40, "Base Mesh");
      o->box(FL_GTK_UP_BOX);
      o->down_box(FL_GTK_DOWN_BOX);
      o->value(1);
      o->selection_color((Fl_Color)1);
      o->callback((Fl_Callback*)cb_Base);
    } // Fl_Light_Button* o
    { Fl_Button* o = new Fl_Button(1010, 180, 85, 40, "Texture");
      o->box(FL_GTK_UP_BOX);
      o->down_box(FL_GTK_DOWN_BOX);
      o->callback((Fl_Callback*)cb_Texture);
    } // Fl_Button* o
    { Fl_Button* o = new Fl_Button(1110, 180, 85, 40, "Environment");
      o->box(FL_GTK_UP_BOX);
      o->down_box(FL_GTK_DOWN_BOX);
      o->callback((Fl_Callback*)cb_Environment);
    } // Fl_Button* o
    { Fl_Light_Button* o = new Fl_Light_Button(1010, 225, 85, 40, "   @+6>");
      o->box(FL_GTK_UP_BOX);
      o->down_box(FL_GTK_DOWN_BOX);
      o->selection_color((Fl_Color)1);
      o->callback((Fl_Callback*)cb_6);
    } // Fl_Light_Button* o
    { Fl_Value_Slider* o = new Fl_Value_Slider(1022, 354, 50, 205, "Displace");
      o->type(4);
      o->box(FL_GTK_DOWN_BOX);
      o->color((Fl_Color)FL_BACKGROUND2_COLOR);
      o->selection_color((Fl_Color)1);
      o->value(0.05);
      o->textsize(14);
      o->callback((Fl_Callback*)cb_Displace);
      o->align(FL_ALIGN_TOP);
    } // Fl_Value_Slider* o
    { Fl_Value_Slider* o = new Fl_Value_Slider(1025, 591, 50, 205, "Width");
      o->type(4);
      o->box(FL_GTK_DOWN_BOX);
      o->color((Fl_Color)FL_BACKGROUND2_COLOR);
      o->selection_color((Fl_Color)1);
      o->maximum(1.5);
      o->step(0.05);
      o->value(0.5);
      o->textsize(14);
      o->callback((Fl_Callback*)cb_Width);
      o->align(FL_ALIGN_TOP);
    } // Fl_Value_Slider* o
    { Fl_Value_Slider* o = new Fl_Value_Slider(1125, 593, 50, 205, "Curvature");
      o->type(4);
      o->box(FL_GTK_DOWN_BOX);
      o->color((Fl_Color)FL_BACKGROUND2_COLOR);
      o->selection_color((Fl_Color)1);
      o->maximum(2);
      o->step(0.05);
      o->value(0.5);
      o->textsize(14);
      o->callback((Fl_Callback*)cb_Curvature);
      o->align(FL_ALIGN_TOP);
    } // Fl_Value_Slider* o
    { Fl_Value_Slider* o = new Fl_Value_Slider(1140, 30, 35, 145, "center");
      o->type(4);
      o->box(FL_GTK_DOWN_BOX);
      o->color((Fl_Color)FL_BACKGROUND2_COLOR);
      o->selection_color((Fl_Color)1);
      o->minimum(0.005);
      o->step(0.05);
      o->value(1);
      o->textsize(14);
      o->callback((Fl_Callback*)cb_center);
      o->align(FL_ALIGN_TOP);
    } // Fl_Value_Slider* o
    { Fl_Group* o = new Fl_Group(1100, 251, 116, 304);
      o->box(FL_EMBOSSED_FRAME);
      o->color((Fl_Color)6);
      { Fl_Button* o = new Fl_Button(1110, 445, 85, 40, "flatten");
        o->box(FL_GTK_UP_BOX);
        o->down_box(FL_GTK_DOWN_BOX);
        o->callback((Fl_Callback*)cb_flatten);
      } // Fl_Button* o
      { Fl_Button* o = new Fl_Button(1110, 495, 85, 40, "write_file");
        o->box(FL_GTK_UP_BOX);
        o->down_box(FL_GTK_DOWN_BOX);
        o->callback((Fl_Callback*)cb_write_file);
      } // Fl_Button* o
      { Fl_Value_Input* o = new Fl_Value_Input(1113, 406, 75, 24, "hole radius");
        o->minimum(0.3);
        o->step(0.05);
        o->value(0.05);
        o->callback((Fl_Callback*)cb_hole);
        o->align(FL_ALIGN_TOP);
      } // Fl_Value_Input* o
      { Fl_Value_Input* o = new Fl_Value_Input(1112, 280, 75, 24, "paper width");
        o->minimum(20);
        o->maximum(128);
        o->step(0.5);
        o->value(16.5);
        o->callback((Fl_Callback*)cb_paper);
        o->align(FL_ALIGN_TOP);
      } // Fl_Value_Input* o
      { Fl_Value_Input* o = new Fl_Value_Input(1112, 322, 75, 24, "paper height");
        o->minimum(10);
        o->maximum(128);
        o->step(0.5);
        o->value(22);
        o->callback((Fl_Callback*)cb_paper1);
        o->align(FL_ALIGN_TOP);
      } // Fl_Value_Input* o
      { Fl_Value_Input* o = new Fl_Value_Input(1115, 361, 75, 24, "piece size");
        o->minimum(5);
        o->maximum(12);
        o->step(0.1);
        o->value(5);
        o->callback((Fl_Callback*)cb_piece);
        o->align(FL_ALIGN_TOP);
      } // Fl_Value_Input* o
      o->end();
    } // Fl_Group* o
    o->end();
  } // Fl_Double_Window* o
  Fl::visual(FL_DOUBLE | FL_INDEX);
  w->show(argc, argv);
  return Fl::run();
}
