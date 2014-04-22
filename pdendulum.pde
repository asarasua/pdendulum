/*
* Alvaro Sarasua
* Double pendulum
* 
*/ 
DoublePendulum pendulum = new DoublePendulum (0.005, 9.8);
final int THETA_1 = 0;
final int OMEGA_1 = 1;
final int THETA_2 = 2;
final int OMEGA_2 = 3;
final int NUM_EQNS = 4;
final int SHADE_SIZE = 250;

boolean ending = false;
boolean phoneTapped = false;
int tapTime = 0;

ArrayList<PVector> shade;

void setup() {
  size(192, 157, P3D); 
  background(0); 
  frameRate(25);

  shade = new ArrayList<PVector>();
  
  //Initialize pendulum 
  pendulum.theta1 = 0.8;
  pendulum.omega1 = 0.0;
  pendulum.m1 = 0.9;
  pendulum.l1 = 0.9;
    
  pendulum.theta2 = 0.5 + random(-0.0001, 0.0001);
  pendulum.omega2 = 0.0;
  pendulum.m2 = 0.1;
  pendulum.l2 = 0.45;
  
  phones.connect();
}

void draw () {
  pendulum.updateTime(millis());
  
  background(0);
  lights();

  float xoff = 0;
  
  pushMatrix();
    //Draw pendulum
    translate(width/2, 0);

    float x1 = pendulum.l1 * sin(pendulum.theta1);
    float y1 = pendulum.l1 * cos(pendulum.theta1);
    float x2 = x1 + pendulum.l2 * sin(pendulum.theta2);
    float y2 = y1 + pendulum.l2 * cos(pendulum.theta2);
    
    x1 *= 100;
    y1 *= 100;
    x2 *= 100;
    y2 *= 100;

    //update shade
    if (millis() > 2000){
      shade.add(new PVector(x2, y2));
    }
    
    if (shade.size() >= SHADE_SIZE) {
      shade.remove(0);
    }

    //draw shade
    for (int i = 0; i < shade.size()-1; ++i) {
      //250, 255, 22 yellow
      stroke(250, 255 * (float)i / shade.size(), 22, map(i, 0, SHADE_SIZE, 10, 255));
      strokeWeight(5 * (float)i / shade.size());
      randomizer = millis()/700 * (float)(shade.size()-i) / shade.size();
      if (phoneTapped){
        randomizer *= 10;
        if (millis() > tapTime + 500){
          phoneTapped = false;
        }
      }
      if (ending){
        randomizer = 1000;
      }
      xoff += 0.05;
      
      //line(shade.get(i).x + random(-randomizer,randomizer), shade.get(i).y + random(-randomizer,randomizer), (float)i / shade.size(), shade.get(i+1).x + random(-randomizer,randomizer), shade.get(i+1).y + random(-randomizer,randomizer), (float)(i+1) / shade.size());
      line(shade.get(i).x + map(noise(xoff), 0, 1, -randomizer,randomizer), shade.get(i).y + map(noise(xoff), 0, 1, -randomizer,randomizer), (float)i / shade.size(), shade.get(i+1).x + map(noise(xoff), 0, 1, -randomizer,randomizer), shade.get(i+1).y + map(noise(xoff), 0, 1, -randomizer,randomizer), (float)(i+1) / shade.size());
    }

    //stroke(130, 82, 1);
    stroke(200, 200, 200);
    strokeWeight(0.7);
    line(0,0, x1, y1);
    line(x1, y1, x2, y2);
    noStroke();

    pushMatrix();  
      translate(x1, y1, 0);
      fill(201, 192, 187);
      sphere(3);
      //add light
    popMatrix();
    pushMatrix();
      translate(x2, y2, 0);
      //ellipse(x1, y1, 10, 10);
      //fill(229, 228, 226);
      fill(131,156,165);
      sphere(5);
    popMatrix();
    //ellipse(x2, y2, 10, 10);

  popMatrix();
}


class DoublePendulum{
  /**
     * Angle of the first pendulum from the vertical (in rad).
     */
    float theta1;

    /**
     * Angular acceleration of the first pendulum (dθ/dt).
     */
    float omega1;

    /**
     * Length of the first pendulum (in m).
     */
    float l1;

    /**
     * Mass of the first pendulum (in kg).
     */
    float m1;

    /**
     * Angle of the second pendulum from the vertical (in rad).
     */
    float theta2;

    /**
     * Angular acceleration of the second pendulum (dθ/dt).
     */
    float omega2;

    /**
     * Length of the second pendulum (in m).
     */
    float l2;

    /**
     * Mass of the second pendulum (in m).
     */
    float m2;

    /**
     * Step size to take when numerically solving the ODE.
     */
    float dt;

    /**
     * Acceleration due to gravity (usually 9.81 ms^-2).
     */
    float g;

    /**
     * Current time for which omega and theta are evaluated for.
     */
    float eTime;
    
    DoublePendulum (float _dt, float _g){
      dt = _dt;
      g = _g;
    }

    void updateTime(int newTime){
      float actualTime = newTime / 1000.0;
      if (actualTime > eTime) {
        update(actualTime);
      }
    }

    void update(float newTime){
      do
      {
          float yin[] = { theta1, omega1, theta2, omega2 };
          float yout[] = new float[NUM_EQNS];

          solveODEs(yin, yout);

          theta1 = yout[THETA_1];
          omega1 = yout[OMEGA_1];
          theta2 = yout[THETA_2];
          omega2 = yout[OMEGA_2];
      } while ((eTime += dt) < newTime);
    }

    void solveODEs(float[] yin, float[] yout){
      float dydx[] = new float[NUM_EQNS];
      float dydxt[] =  new float[NUM_EQNS];
      float yt[] = new float[NUM_EQNS];
      float k1[] = new float[NUM_EQNS];
      float k2[] = new float[NUM_EQNS];
      float k3[] = new float[NUM_EQNS];
      float k4[] = new float[NUM_EQNS];

      // First step
      derivs(yin, dydx);
      for (int i = 0; i < NUM_EQNS; ++i)
      {
          k1[i] = dt * dydx[i];
          yt[i] = yin[i] + 0.5 * k1[i];
      }

      // Second step
      derivs(yt, dydxt);
      for (int i = 0; i < NUM_EQNS; ++i)
      {
          k2[i] = dt * dydxt[i];
          yt[i] = yin[i] + 0.5 * k2[i];
      }

      // Third step
      derivs(yt, dydxt);
      for (int i = 0; i < NUM_EQNS; ++i)
      {
          k3[i] = dt * dydxt[i];
          yt[i] = yin[i] + k3[i];
      }

      // Fourth step
      derivs(yt, dydxt);
      for (int i = 0; i < NUM_EQNS; ++i)
      {
          k4[i] = dt * dydxt[i];
          yout[i] = yin[i] + k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0;
      }
    }

  void derivs(float[] yin, float[] dydx){
    // Delta is θ2 - θ1
    float delta = yin[THETA_2] - yin[THETA_1];

    // `Big-M' is the total mass of the system, m1 + m2;
    float M = m1 + m2;

    // Denominator expression for ω1
    float den = M*l1 - m2*l1*cos(delta)*cos(delta);

    // dθ/dt = ω, by definition
    dydx[THETA_1] = yin[OMEGA_1];

    // Compute ω1
    dydx[OMEGA_1] = (m2*l1*yin[OMEGA_1]*yin[OMEGA_1]*sin(delta)*cos(delta)
                  + m2*g*sin(yin[THETA_2])*cos(delta)
                  + m2*l2*yin[OMEGA_2]*yin[OMEGA_2]*sin(delta)
                  - M*g*sin(yin[THETA_1])) / den;

    // Again, dθ/dt = ω for θ2 as well
    dydx[THETA_2] = yin[OMEGA_2];

    // Multiply den by the length ratio of the two bobs
    den *= l2 / l1;

    // Compute ω2
    dydx[OMEGA_2] = (-m2*l2*yin[OMEGA_2]*yin[OMEGA_2]*sin(delta)*cos(delta)
                  + M*g*sin(yin[THETA_1])*cos(delta)
                  - M*l1*yin[OMEGA_1]*yin[OMEGA_1]*sin(delta)
                  - M*g*sin(yin[THETA_2])) / den;
  }
}

phones.onClick(function(id) {
    phoneTapped = true;
    tapTime = millis();
});

sketch.onAboutToStop(function() { 
  ending = true;
});
