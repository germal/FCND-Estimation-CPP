#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"
#include <time.h>

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
  //counter = 0;
  //lastt = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config + ".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  kiPosZ = config->Get(_config + ".kiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // Access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the value kappa [m]:
  // torque (Nm) produced by motor per N of thrust generated.

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float idealThrust = collThrustCmd / 4.f;
  float d = 4.f * L * M_SQRT1_2; // adjust L to be l
  float p = momentCmd.x / d;
  float q = momentCmd.y / d;
  float r = momentCmd.z / (kappa*4.f);

  cmd.desiredThrustsN[0] = idealThrust + p + q - r; // front left  F1
  cmd.desiredThrustsN[1] = idealThrust - p + q + r; // front right F2
  cmd.desiredThrustsN[2] = idealThrust + p - q + r; // rear left   F4
  cmd.desiredThrustsN[3] = idealThrust - p - q - r; // rear right  F3

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  V3F I = V3F(Ixx, Iyy, Izz);
  // Account for the moment of inertia of the drone
  momentCmd = I * kpPQR * (pqrCmd - pqr);
  // Gate the moment if necessary
  float maxTorque = 6.28f;
  if (momentCmd.mag() > maxTorque)
    momentCmd = momentCmd * maxTorque / momentCmd.mag();

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0.f;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float err_z = posZCmd - posZ; // Proportion p component
  integratedAltitudeError += err_z * dt;
  // Experiment: prevent integral wind-up
  //if (abs(integratedAltitudeError) > 0.25)
  //  integratedAltitudeError = 0;
  float newAcc = kpPosZ * err_z + kpVelZ * (velZCmd - velZ) // add deriv d component
    + kiPosZ * integratedAltitudeError + accelZCmd;         // and i component
  // Gate the acceleration if necessary
  newAcc = CONSTRAIN(newAcc, -maxAscentRate / dt, maxDescentRate / dt);
  float acc = (newAcc - float(CONST_GRAVITY)) / R(2,2);
  thrust = -acc * mass; // convert to thrust force [N]

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return thrust;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  // Account for drone mass (factor out for acceleration)
  float acc = collThrustCmd / mass;
  // Gate the tilt angle (approximately)
  float bx_targ = -CONSTRAIN(accelCmd.x / acc, -maxTiltAngle, maxTiltAngle);
  float by_targ = -CONSTRAIN(accelCmd.y / acc, -maxTiltAngle, maxTiltAngle);
  // Calculate angular velocity
  float bxd = kpBank * (bx_targ - R(0,2));
  float byd = kpBank * (by_targ - R(1,2));
  pqrCmd.x = (R(1,0) * bxd - R(1,1) * byd) / R(2,2);
  pqrCmd.y = (R(0,0) * bxd - R(0,1) * byd) / R(2,2);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to constrain float foo to range [0,b]
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd = 0.f;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  // Constrain to valid yaw range of -pi to pi
  float errYaw = yawCmd - yaw;
  if (errYaw > float(M_PI))
    errYaw = errYaw - 2.f * float(M_PI);
  else if (errYaw < -float(M_PI))
    errYaw = errYaw + 2.f * float(M_PI);
  yawRateCmd = kpYaw * errYaw;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmd)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmd: desired acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use fmodf(foo,b) to constrain float foo to range [0,b]
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you cap the horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmd.z = 0.f;
  velCmd.z = 0.f;
  posCmd.z = pos.z;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  V3F kpP(kpPosXY, kpPosXY, 0.f); // do not apply gain to z-component
  V3F kdV(kpVelXY, kpVelXY, 0.f);
  // Not found to be needed here, but: constrain velocity target
  //if (velCmd.mag() > maxSpeedXY)
  //  velCmd = velCmd * maxSpeedXY / velCmd.mag();

  // Compute acceleration
  accelCmd += kpP * (posCmd - pos) + kdV * (velCmd - vel);

  // Constrain computed acceleration
  if (accelCmd.mag() > maxAccelXY)
    accelCmd = accelCmd * maxAccelXY / accelCmd.mag();

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
/*int t = time(NULL);
  counter += 1;
  if (t != lastt) {
    cout << counter << endl;
    counter = 0;
  }
  lastt = t; */
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
