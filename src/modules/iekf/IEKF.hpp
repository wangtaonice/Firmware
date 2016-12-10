#include "ros/ros.hpp"
#include "matrix/math.hpp"

using namespace matrix;

/**
 * Note that structs are used instead of enums
 * to allow use in arrays without casting
 * and to keep size small
 */

/**
 * State enum
 */
struct X {
	static const uint8_t q_nb_0 = 0;
	static const uint8_t q_nb_1 = 1;
	static const uint8_t q_nb_2 = 2;
	static const uint8_t q_nb_3 = 3;
	static const uint8_t vel_n = 4;
	static const uint8_t vel_e = 5;
	static const uint8_t vel_d = 6;
	static const uint8_t gyro_bias_bx = 7;
	static const uint8_t gyro_bias_by = 8;
	static const uint8_t gyro_bias_bz = 9;
	static const uint8_t accel_scale = 10;
	static const uint8_t pos_n = 11;
	static const uint8_t pos_e = 12;
	static const uint8_t pos_d = 13;
	static const uint8_t terrain_alt = 14;
	static const uint8_t baro_bias = 15;
	static const uint8_t n = 16;
};

/**
 * Error state enum
 * used for linearization
 *
 * Note gyro bias in navigation frame
 */
struct Xe {
	static const uint8_t rot_bx = 0;
	static const uint8_t rot_by = 1;
	static const uint8_t rot_bz = 2;
	static const uint8_t vel_n = 3;
	static const uint8_t vel_e = 4;
	static const uint8_t vel_d = 5;
	static const uint8_t gyro_bias_n = 6;
	static const uint8_t gyro_bias_e = 7;
	static const uint8_t gyro_bias_d = 8;
	static const uint8_t accel_scale = 9;
	static const uint8_t pos_n = 10;
	static const uint8_t pos_e = 11;
	static const uint8_t pos_d = 12;
	static const uint8_t terrain_alt = 13;
	static const uint8_t baro_bias = 14;
	static const uint8_t n = 15;
};

/**
 * Input enum
 */
struct U {
	static const uint8_t omega_nb_bx = 0;
	static const uint8_t omega_nb_by = 1;
	static const uint8_t omega_nb_bz = 2;
	static const uint8_t accel_bx = 3;
	static const uint8_t accel_by = 4;
	static const uint8_t accel_bz = 5;
	static const uint8_t n = 6;
};

/**
 * Accel measurement enum
 */
struct Y_accel {
	static const uint8_t accel_bx = 0;
	static const uint8_t accel_by = 1;
	static const uint8_t accel_bz = 2;
	static const uint8_t n = 3;
};

/**
 * GPS measurement
 */
struct Y_gps {
	static const uint8_t pos_n = 0;
	static const uint8_t pos_e = 1;
	static const uint8_t pos_d = 2;
	static const uint8_t vel_n = 3;
	static const uint8_t vel_e = 4;
	static const uint8_t vel_d = 5;
	static const uint8_t n = 6;
};

/**
 * Baro measurement
 */
struct Y_baro {
	static const uint8_t asl = 0;
	static const uint8_t n = 1;
};

/**
 * Magnetometer measurement
 *
 * The filter treats the error
 * in the navigation frame
 * (north, east, down) even though the
 * field is measured in the body
 * frame.
 */
struct Y_mag {
	static const uint8_t mag_n = 0;
	static const uint8_t mag_e = 1;
	static const uint8_t mag_d = 2;
	static const uint8_t n = 3;
};

static const float BETA_TABLE[] = {
	0,
	8.82050518214,
	12.094592431,
	13.9876612368,
	16.0875642296,
	17.8797700658,
	19.6465647819,
	21.3802576894,
	23.0806434845,
	24.6673803845,
	26.1487953661,
	27.6350821245,
	29.6565383703,
	31.2211113844,
	32.7673547211,
	34.2967756977,
	35.6906782236,
	37.0724753352,
	38.4549693067,
	39.836592699,
};

/**
 * Main class for invariant extended kalman filter
 *
 * inspired by: https://hal.archives-ouvertes.fr/hal-00494342/document
 *
 * Also see python directory for simulation and python version.
 */
class IEKF
{
public:
	IEKF();
	Vector<float, X::n> dynamics(const Vector<float, X::n> &x, const Vector<float, U::n> &u);
	bool ok() { return _nh.ok(); }
	void callback_gyro(const sensor_gyro_s *msg);
	void callback_accel(const sensor_accel_s *msg);
	void callback_mag(const sensor_mag_s *msg);
	void callback_baro(const sensor_baro_s *msg);
	void callback_gps(const vehicle_gps_position_s *msg);
	void predict(float dt);

	/**
	 * Measurement correction
	 *
	 * @param r : measurement residual, r =  y - g(x)
	 * @param H : measurement matrix y = g(x, u), H = dg(x, u)/dx
	 * @param R : measurement covariance matrix R = E[r^Tr]
	 * @return : true, no fault (correction applied), false (fault detect, did not correct)
	 **/
	template<size_t n_y>
	bool correct(Vector<float, n_y> &r, Matrix<float, n_y, Xe::n> &H,
		     Matrix<float, n_y, n_y> &R)
	{
		SquareMatrix<float, n_y> S_I = SquareMatrix<float, n_y>(H * _P * H.T() + R).I();
		Matrix<float, Xe::n, n_y> K = _P * H.T() * S_I;

		float beta = (r.T() * (S_I * r))(0, 0);
		bool fault = beta > BETA_TABLE[n_y];
		//ROS_INFO("beta: %10.4f, thresh: %10.4f\n", double(beta), double(BETA_TABLE[n_y]));

		if (!fault) {
			// TODO apply all correction
			Vector<float, Xe::n> d_xe = K * r;
			Quaternion<float> q_nb(_x(X::q_nb_0), _x(X::q_nb_1), _x(X::q_nb_2), _x(X::q_nb_3));
			Quaternion<float> d_q_nb = Quaternion<float>(0,
						   d_xe(Xe::rot_bx), d_xe(Xe::rot_by), d_xe(Xe::rot_bz)) * q_nb;
			//ROS_INFO("d_q_nb");
			//d_q_nb.print();
			Vector3<float> d_gyro_bias_b = q_nb.conjugate_inversed(
							       Vector3<float>(d_xe(Xe::gyro_bias_n),
									       d_xe(Xe::gyro_bias_e),
									       d_xe(Xe::gyro_bias_d)));

			// linear term correction is the same
			// as the error correction
			_x(X::q_nb_0) += d_q_nb(0);
			_x(X::q_nb_1) += d_q_nb(1);
			_x(X::q_nb_2) += d_q_nb(2);
			_x(X::q_nb_3) += d_q_nb(3);
			_x(X::gyro_bias_bx) += d_gyro_bias_b(0);
			_x(X::gyro_bias_by) += d_gyro_bias_b(1);
			_x(X::gyro_bias_bz) += d_gyro_bias_b(2);
			_x(X::vel_n) += d_xe(Xe::vel_n);
			_x(X::vel_e) += d_xe(Xe::vel_e);
			_x(X::vel_d) += d_xe(Xe::vel_d);
			_x(X::pos_n) += d_xe(Xe::pos_n);
			_x(X::pos_e) += d_xe(Xe::pos_e);
			_x(X::pos_d) += d_xe(Xe::pos_d);
			_x(X::terrain_alt) += d_xe(Xe::terrain_alt);
			_x(X::baro_bias) += d_xe(Xe::baro_bias);
			_P -= K * H * _P;
		}

		return fault;
	};

private:
	ros::NodeHandle _nh;
	ros::Subscriber _sub_gyro;
	ros::Subscriber _sub_accel;
	ros::Subscriber _sub_mag;
	ros::Subscriber _sub_baro;
	ros::Subscriber _sub_attitude;
	ros::Publisher _pub_attitude;
	ros::Publisher _pub_local_position;
	ros::Publisher _pub_global_position;
	Vector<float, X::n> _x; // state
	Matrix<float, Xe::n, Xe::n> _P; // covariance
	Vector<float, U::n> _u;
	Vector3<float> _g_n;
	Vector3<float> _B_n;
};
