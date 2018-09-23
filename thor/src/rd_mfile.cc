/* Author:       Johan Bengtsson

   Definitions:  Procedure for reading a flat file.              */

/* Procedure to read machine file.

   Type codes:

     marker     -1
     drift	 0
     multipole   1
     cavity      2
     thin kick   3
     wiggler     4
     kick map    6

   Integration methods:

     fourth order symplectic integrator   4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps
     apertures: xmin, xmax, ymin, ymax

   The following lines follows depending on element type.

     type

     drift:	 L

     multipole:  hor., ver. displ., roll angle (design), roll angle (error)
                 L, 1/rho, entrance angle, exit angle
                 apertures[4]
		 no of nonzero multipole coeff., n design
		 n, b , a
		     n   n
		     .
		     .
		     .

     wiggler:    L [m], lambda [m], BoBrho [m^-1], k_x [m]

     cavity:	 cavity voltage/beam energy [eV], omega/c, beam energy [eV],
                 phi

     thin kick:	 hor., ver. displacement, roll angle (total)
		 no of nonzero multipole coeff.
		 n, b , a
		     n   n
		     .
		     .
		     .

     kick_map:   scale order <file name>

*/


const int line_max = 200;


template<typename T>
void clr_mpole(mpole_type<T> *mpole)
{
  int  i;

  mpole->method = 0; mpole->n_step = 0;
  mpole->dx_sys[X_] = 0.0; mpole->dx_sys[Y_] = 0.0;
  mpole->dx_rms[X_] = 0.0; mpole->dx_rms[Y_] = 0.0;
  mpole->dx_rnd[X_] = 0.0; mpole->dx_rnd[Y_] = 0.0;
  mpole->droll_par = 0.0; mpole->droll_sys = 0.0; mpole->droll_rms = 0.0;
  mpole->droll_rnd = 0.0;
  for (i = 0; i < mpole_max; i++) {
    mpole->an_par[i] = 0.0; mpole->an_sys[i] = 0.0; mpole->an_rms[i] = 0.0;
    mpole->an_rnd[i] = 0.0; mpole->an[i] = 0.0;
    mpole->bn_par[i] = 0.0; mpole->bn_sys[i] = 0.0; mpole->bn_rms[i] = 0.0;
    mpole->bn_rnd[i] = 0.0; mpole->bn[i] = 0.0;
  }
  mpole->order = 0; mpole->n_design = 0; mpole->edge1 = 0.0;
  mpole->edge2 = 0.0; mpole->gap = 0.0; mpole->h_bend = 0.0;
}

template<typename T>
void clr_elem(elem_type<T> &elem)
{
  elem.L = 0.0; elem.S = 0.0; elem.Fnum = 0; elem.Knum = 0;
  elem.dx[X_] = 0.0; elem.dx[Y_] = 0.0;
  elem.droll[X_] = 0.0; elem.droll[Y_] = 0.0; elem.droll_par = 0.0;
  elem.c0 = 0.0; elem.c1 = 0.0; elem.s1 = 0.0;
  elem.kind = Undef;
  elem.max_ampl[X_][0] = 0.0; elem.max_ampl[X_][0] = 0.0;
  elem.max_ampl[Y_][0] = 0.0; elem.max_ampl[Y_][0] = 0.0;
}

void clr_Family(void)
{
  int  i, j;

  for (i = 0; i < max_Family; i++) {
    Families[i].n_Kids = 0;
    for (j = 0; j < max_Kids; j++)
      Families[i].Kids[j] = 0;
  }
}

template<typename T>
void rd_mfile(const char file_name[], elem_type<T> elem[])
{
  const int n_ps = 6;

  long int       ind;
  int            j, k, order, n_mpole, Fnum, Knum, method, n_step;
  double         drerror, bn, an, L, val[n_ps];
  char           line[line_max];
  ss_vect<tps>   Id;
  std::ifstream  inf;

  bool  prt = false;

  std::cout << std::endl;
  std::cout << "reading machine file: " << file_name <<std:: endl;;

  file_rd(inf, file_name);

  clr_Family(); n_elem = 0;
  while (inf.getline(line, line_max)) {
    n_elem++;
    if (n_elem <= max_elem) {
      sscanf(line, "%*s %*d %*d %ld", &ind);
      clr_elem(elem[ind]);

      sscanf(line, "%s %d %d", elem[ind].Name, &Fnum, &Knum);
      k = strlen(elem[ind].Name);
      if (k >= name_length) {
	printf("rd_mfile: name_length exceeded %d (%d) %s\n",
	       k, name_length, elem[ind].Name);
	exit(1);
      }
      elem[ind].Fnum = Fnum; elem[ind].Knum = Knum;
      if ((0 < Fnum) && (Fnum <= max_Family)) {
	strcpy(Families[Fnum-1].Name, elem[ind].Name);
	if ((0 < Knum) && (Knum <= max_Kids)) {
	  Families[Fnum-1].Kids[Knum-1] = ind;
	  Families[Fnum-1].n_Kids = max(Knum, Families[Fnum-1].n_Kids);
	} else if (0 < Knum) {
	  printf("rd_mfile: max_Kids exceeded %d (%d)\n", Knum, max_Kids);
	  exit(1);
	}
      } else if (0 < Fnum) {
	printf("rd_mfile: max_Family exceeded %d (%d)\n", Fnum, max_Family);
	exit(1);
      }

      inf.getline(line, line_max);
      sscanf(line, "%d %d %d", &elem[ind].kind, &method, &n_step);
      if (prt)
	std::cout << elem[ind].Name << " " << elem[ind].kind << std::endl;

      inf.getline(line, line_max); 
      sscanf(line, "%lf %lf %lf %lf",
	     &elem[ind].max_ampl[X_][0], &elem[ind].max_ampl[X_][1],
	     &elem[ind].max_ampl[Y_][0], &elem[ind].max_ampl[Y_][1]);
      if (prt) std::cout << std::scientific << std::setprecision(3)
	     << std::setw(11) << elem[ind].max_ampl[X_][0]
	     << std::setw(11) << elem[ind].max_ampl[X_][0]
	     << std::setw(11) << elem[ind].max_ampl[Y_][1]
	     << std::setw(11) << elem[ind].max_ampl[Y_][1] << std::endl;

      switch (elem[ind].kind) {
      case Marker:
	L = 0.0;
	break ;
      case Drift:
	// Note, L is polymorphic
	inf.getline(line, line_max); sscanf(line, "%lf", &L);
	break;
      case Mpole:
	elem[ind].mpole = new mpole_type<T>;

	clr_mpole(elem[ind].mpole);

	elem[ind].mpole->method = method;
	elem[ind].mpole->n_step = n_step;

	inf.getline(line, line_max);
	sscanf(line, "%lf %lf %lf %lf",
	       &elem[ind].dx[X_], &elem[ind].dx[Y_],
	       &elem[ind].mpole->droll_par, &drerror);
	elem[ind].droll[X_] =
	  cos(dtor(elem[ind].mpole->droll_par+drerror));
	elem[ind].droll[Y_] =
	  sin(dtor(elem[ind].mpole->droll_par+drerror));

	inf.getline(line, line_max);
	sscanf(line, "%lf %lf %lf %lf %lf",
	       &L, &elem[ind].mpole->h_bend,
	       &elem[ind].mpole->edge1, &elem[ind].mpole->edge2,
	       &elem[ind].mpole->gap);
	elem[ind].c0 = sin(L*elem[ind].mpole->h_bend/2.0);
	elem[ind].c1 = cos(dtor(elem[ind].mpole->droll_par))
	  *elem[ind].c0;
	elem[ind].s1 = sin(dtor(elem[ind].mpole->droll_par))
	  *elem[ind].c0;
	elem[ind].mpole->droll_rms = drerror - elem[ind].mpole->droll_par;
	elem[ind].mpole->droll_rnd = 1e0;

	inf.getline(line, line_max);
	sscanf(line, "%d %d", &n_mpole, &elem[ind].mpole->n_design);
	for (k = 1; k <= n_mpole; k++) {
	  inf.getline(line, line_max);
	  sscanf(line, "%d %lf %lf", &order, &bn, &an);
	  if (order <= mpole_max) {
	    elem[ind].mpole->bn[order-1] = bn;
	    elem[ind].mpole->an[order-1] = an;
	    elem[ind].mpole->order = order;
	  } else {
	    std::cout << "rd_mfile: mpole_max exceeded " << order
		 << "(" << mpole_max << ")" << std::endl;
	    exit(1);
	  }
	}
	break;
      case Cavity:
	elem[ind].cavity = new cavity_type;

	inf.getline(line, line_max);
	sscanf(line, "%lf %lf %d %lf %lf",
	       &elem[ind].cavity->V_rf, &elem[ind].cavity->f_rf,
	       &elem[ind].cavity->h_rf, &E0, &elem[ind].cavity->phi_rf);
	elem[ind].cavity->V_rf = elem[ind].cavity->V_rf*E0;
	elem[ind].cavity->f_rf = elem[ind].cavity->f_rf*clight/(2.0*pi);
	E0 = E0/1e9; L = 0.0;
	break;
      case Thinkick:
	elem[ind].kind = Mpole;
	elem[ind].mpole = new mpole_type<T>;

	clr_mpole(elem[ind].mpole);

	elem[ind].mpole->method = method;
	elem[ind].mpole->n_step = n_step;

	inf.getline(line, line_max);
	sscanf(line, "%lf %lf %lf", &elem[ind].dx[X_],
	       &elem[ind].dx[Y_], &drerror);

	L = 0.0;
	elem[ind].droll[X_] = cos(dtor(drerror));
	elem[ind].droll[Y_] = sin(dtor(drerror));
	elem[ind].c0 = sin(L*elem[ind].mpole->h_bend/2.0);
	elem[ind].c1 =
	  cos(dtor(elem[ind].mpole->droll_par))*elem[ind].c0;
	elem[ind].s1 =
	  sin(dtor(elem[ind].mpole->droll_par))*elem[ind].c0;
	elem[ind].mpole->droll_rms = drerror; elem[ind].mpole->droll_rnd = 1e0;

	inf.getline(line, line_max);
	sscanf(line, "%d %d", &n_mpole, &elem[ind].mpole->n_design);
	for (k = 1; k <= n_mpole; k++) {
	  inf.getline(line, line_max);
	  sscanf(line, "%d %lf %lf", &order, &bn, &an);
	  if (order <= mpole_max) {
	    elem[ind].mpole->bn[order-1] = bn;
	    elem[ind].mpole->an[order-1] = an;
	    elem[ind].mpole->order = order;
	  } else {
	    std::cout << "rdmfile: mpole_max exceeded " << order
		 << "(" << mpole_max << ")" << std::endl;
	    exit(0);
	  }
	}
	break;
      case Wiggler:
	elem[ind].wiggler = new wiggler_type;
	elem[ind].wiggler->method = method;
	elem[ind].wiggler->n_step = n_step;

	inf.getline(line, line_max);
	sscanf(line, "%lf %lf", &L, &elem[ind].wiggler->lambda);

	inf.getline(line, line_max);
	sscanf(line, "%d", &elem[ind].wiggler->n_harm);
	for (k = 0; k < elem[ind].wiggler->n_harm; k++) {
	  inf.getline(line, line_max);
	  sscanf(line, "%d %lf %lf %lf %lf %lf",
		 &elem[ind].wiggler->harm[k],
		 &elem[ind].wiggler->kxV[k],
		 &elem[ind].wiggler->BoBrhoV[k],
		 &elem[ind].wiggler->kxH[k],
		 &elem[ind].wiggler->BoBrhoH[k],
		 &elem[ind].wiggler->phi[k]);
	}
	break;
      case Kick_map:
	elem[ind].kick_map = new kick_map_type;
	elem[ind].kick_map->method = method;
	elem[ind].kick_map->n_step = n_step;

	inf.getline(line, line_max);
	sscanf(line, "%lf %d %s", &elem[ind].kick_map->scl,
	       &elem[ind].kick_map->order, elem[ind].kick_map->file_name);

	if (elem[ind].kick_map->order == 1)
	  Read_IDfile(elem[ind].kick_map->file_name, L,
		      elem[ind].kick_map->nx, elem[ind].kick_map->nz,
		      elem[ind].kick_map->tabx, elem[ind].kick_map->tabz,
		      elem[ind].kick_map->thetax1,
		      elem[ind].kick_map->thetaz1);
	else
	  Read_IDfile(elem[ind].kick_map->file_name, L,
		      elem[ind].kick_map->nx, elem[ind].kick_map->nz,
		      elem[ind].kick_map->tabx, elem[ind].kick_map->tabz,
		      elem[ind].kick_map->thetax, elem[ind].kick_map->thetaz);

	if (elem[ind].kick_map->method == 2) {
	  elem[ind].kick_map->tx = dmatrix(1, elem[ind].kick_map->nz, 1,
					   elem[ind].kick_map->nx);
	  elem[ind].kick_map->tz = dmatrix(1, elem[ind].kick_map->nz, 1,
					   elem[ind].kick_map->nx);
	  elem[ind].kick_map->tab1 = (double *)malloc((elem[ind].kick_map->nx)
						      *sizeof(double));
	  elem[ind].kick_map->tab2 = (double *)malloc((elem[ind].kick_map->nz)
						      *sizeof(double));
	  elem[ind].kick_map->f2x = dmatrix(1, elem[ind].kick_map->nz, 1,
					    elem[ind].kick_map->nx);
	  elem[ind].kick_map->f2z = dmatrix(1, elem[ind].kick_map->nz, 1,
					    elem[ind].kick_map->nx);
	  Matrices4Spline(elem[ind].kick_map);

	  // free(tab1); free(tab2);
	  // free_matrix(tx, 1, nz, 1, nx); free_matrix(tz, 1, nz, 1, nx);
	  // free_matrix(f2x, 1, nz, 1, nx); free_matrix(f2z, 1, nz, 1, nx);
      }
      break;
      case Map_:
	elem[ind].map = new map_type;
	Id.identity(); elem[ind].map->M.zero();
	for (j = 0; j < n_ps; j++) {
	  inf.getline(line, line_max);
	  if (prt) printf("%s\n", line);
	  sscanf(line, "%lf %lf %lf %lf %lf %lf",
		 &val[0], &val[1], &val[2], &val[3], &val[4], &val[5]);
	  for (k = 0; k < n_ps; k++)
	    elem[ind].map->M[j] += val[k]*Id[k];
	}
	break;
      default:
	std::cout << "rd_mfile: undefined element " << elem[ind].Name
	     << std::setw(4) << ind << " " << std::setw(2) << elem[ind].kind
	     << std::endl;
	exit(1);
	break;
      }

      elem[ind].L = L;

      if (ind == 0)
	elem[ind].S = 0.0;
      else
	elem[ind].S = elem[ind-1].S + L;
    } else {
      std::cout << "rd_mfile: max_elem exceeded " << n_elem
	   << "(" << max_elem << ")" << std::endl;
      exit(1);
    }
  }

  std::cout << std::endl;
  std::cout << "no of elements:    " << std::setw(12) << n_elem << std::endl;
  std::cout << "circumference [m]: "
       << std::fixed << std::setw(12) << std::setprecision(6)
	    << elem[n_elem-1].S << std::endl;
  std::cout << std::fixed << "Energy [GeV]:        "
       << std::setw(10) << std::setprecision(2) << E0 << std::endl;

  inf.close();
}
