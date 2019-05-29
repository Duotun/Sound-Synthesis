#include "stvk.h"


void  StVKMaterial::first_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const
{
	second_piola_kirchhoff(F, out);
	out = F * out;    //P=F*S, S- Second piola_kirchhoff
}
//Saint Venant-Kirchhoff model for Hyperelastic Materials
void StVKMaterial::second_piola_kirchhoff(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const
{
	Eigen::Matrix3d C, E;
	right_cauchy_green_deformation_tensor(F, C);
	green_strain_tensor(C, E);

	out = Eigen::Matrix3d::Identity()*E.trace() + 2 * mu_*E;
	//out =\lambda*tr(E)*I+2*\MU*E

}
// cauchy=1/J *F*S*F^T, S- second piolar stress tensor, J is the determinant of F
void StVKMaterial::cauchy_stress_tensor(const Eigen::MatrixX3d &F, Eigen::MatrixX3d &out) const
{
	second_piola_kirchhoff(F, out);
	out = F * out*F.transpose();
	out /= F.determinant();

}
double StVKMaterial::strain_energy_density(const Eigen::MatrixX3d &F) const
{
	Eigen::Matrix3d C, E;
	right_cauchy_green_deformation_tensor(F, C);
	green_strain_tensor(C, E);

	//\psi = \mu *tr(E^T*E)+0.5*\lambda*tr(E)^2;
	Eigen::Matrix3d colon = E.transpose()*E;
	const double tr = E.trace()*E.trace();  //tr(E)^2
	return mu_ * colon.trace() + 0.5*lambda_*tr;
}

//linear tetrahedron element - stiff matrix try
void StVKMaterial::stiffness_matrix(const double density, const std::map<int, double> &mass, const std::vector<Eigen::Vector4i> &tetra, const std::vector<Eigen::Vector3d> &positions, Eigen::SparseMatrix<double, Eigen::RowMajor> &stiffmatrix, int pointsize) const
{
	//out should be a sparse matrix
	//positions are for volume, tetra is for index
	//first get volume of each tetrahedron

	std::map<int, double> volume;  //get volume of every tetrahdron from mass
	std::vector<Eigen::Triplet<double>> coeffcients;

	for (int i=0;i<mass.size();i++)
	{
		volume[i] = mass.at(i) / density;
	}

	//a,b,c,d for every tetrahedron 
	for (int cnt = 0; cnt < tetra.size(); cnt++)
	{
		std::vector<double>a, b, c, d;   //permutation
		std::vector<Eigen::Vector3d> pos;
		pos.push_back(positions[tetra[cnt](0)]);
		pos.push_back(positions[tetra[cnt](1)]);
		pos.push_back(positions[tetra[cnt](2)]);    //get four positions
		pos.push_back(positions[tetra[cnt](3)]);

		
		//std::cout << std::endl;

		Eigen::MatrixXd c_triangleview(6, 6); c_triangleview.setZero();
		Eigen::MatrixXd c_matrix(6, 6);   //symmetric  
		c_matrix.setZero(); c_triangleview.setZero();

		//iso right now here
		c_matrix(0, 0) = c_matrix(1, 1) = c_matrix(2, 2) = lambda_ + 2 * mu_;
		c_triangleview(0, 1) = c_triangleview(0, 2) = c_triangleview(1, 2) = lambda_;
		c_matrix(3, 3) = c_matrix(4, 4) = c_matrix(5, 5) = mu_;

		c_matrix += c_triangleview + c_triangleview.transpose();

		//std::cout << "mu: " << c_matrix(5, 5) << std::endl;
		for (int i = 0; i < 4; i++)
		{
			int j, k, l;
			j = (i + 1) % 4; k = (j + 1) % 4; l = (k + 1) % 4;    //index

			//std::cout << "I J K: " << i << j << k<<l;
			Eigen::Matrix3d a_d, b_d, c_d, d_d;   //for determinant,
			a_d.row(0) << pos[j].transpose();  b_d.row(0) << pos[j].transpose();  c_d.row(0) << pos[j].transpose(); d_d.row(0) << pos[j].transpose();
			a_d.row(1) << pos[k].transpose();  b_d.row(1) << pos[k].transpose();  c_d.row(1) << pos[k].transpose(); d_d.row(1) << pos[k].transpose();
			a_d.row(2) << pos[l].transpose();  b_d.row(2) << pos[l].transpose(); c_d.row(2) << pos[l].transpose(); d_d.row(2) << pos[l].transpose();

			
			std::cout << a_d.row(0) << std::endl;
			// 1 -> some column
			Eigen::Vector3d t; t.setIdentity();
			b_d.col(0) << t;  c_d.col(1) << t; d_d.col(2) << t;

			//a_i, b_i, c_i, d_i
			a.push_back(a_d.determinant());
			b.push_back(-b_d.determinant());
			c.push_back(-c_d.determinant());
			d.push_back(-d_d.determinant());
			//std::cout << a_d << std::endl << std::endl;
		}

		//build B matrix  12*6, Ni, parameters ->B input directly
		Eigen::MatrixXd B(6, 12);
		B.setZero();
		//first three
		B(0, 0) = b[0]; B(0, 3) = b[1]; B(0, 6) = b[2]; B(0, 9) = b[3];
		B(1, 0) = c[0]; B(1, 3) = c[1]; B(1, 6) = c[2]; B(1, 9) = c[3];
		B(2, 0) = d[0]; B(2, 3) = d[1]; B(2, 6) = d[2]; B(2, 9) = d[3];

		//last three
		B(3, 0) = c[0]; B(3, 1) = b[0]; B(3, 3) = c[1]; B(3, 4) = b[1]; B(3, 6) = b[2]; B(3, 7) = c[2]; B(3, 9) = c[3]; B(3, 10) = c[3];
		B(4, 1) = d[0]; B(4, 2) = c[0]; B(4, 4) = d[1]; B(4, 5) = c[1]; B(4, 7) = d[2]; B(4, 8) = c[2]; B(4, 10) = d[3]; B(4, 11) = c[3];
		B(5, 0) = d[0]; B(5, 2) = b[0]; B(5, 3) = d[1]; B(5, 5) = b[1]; B(5, 6) = d[2]; B(5, 8) = b[2]; B(5, 9) = d[3]; B(5, 11) = b[3];

		B = B * 1 / 6 / volume[cnt];
		//build k_e = V_e *B^T*C*B
		Eigen::MatrixXd K_e(12, 12);
		K_e = volume[cnt] * B.transpose()*c_matrix*B;

		//triplet building acoording to the positions in global stiff matrix
		for (int nid = 0; nid < 4; nid++)
		{
			int rowid = tetra[cnt][nid];
			for (int nid2 = 0; nid2 < 4; ++nid2)
			{
				int colid = tetra[cnt][nid2] * 3;  //draw 12 * 12
				for (int dir = 0; dir<3; dir++)
					for (int dir2 = 0; dir2<3; dir2++)
					{
						int srow = nid * 3 + dir;
						int scol = nid2 * 3 + dir2;   //full matrix not only upper triangle view
						coeffcients.push_back(Eigen::Triplet<double>(rowid + dir, colid + dir2, K_e(srow, scol)));
					}
			}
		}

	}
	// build from triplet
	stiffmatrix.resize(pointsize * 3, pointsize * 3);
	stiffmatrix.setFromTriplets(coeffcients.begin(), coeffcients.end());




}
