/* Author:  Johan Bengtsson */

/* To generate a lattice flat file.

   Type codes:

     marker     -1
     drift	 0
     multipole   1
     cavity      2
     thin kick   3
     wiggler     4

   Integration methods:

     not applicable                    0
     2nd order symplectic integrator   2
     4th order symplectic integrator   4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps
     apertures: xmin, xmax, ymin, ymax

   The following lines follows depending on element type.

     type

     drift:	 L

     multipole:  hor., ver. displ., roll angle (design),
                                          roll angle (error)
                 L, 1/rho, entrance angle, exit angle
                 apertures[4]
		 no of nonzero multipole coeff., n design
		 n, b , a
		     n   n
		    ...

     wiggler:    L [m], lambda [m]
                 no of harmonics
                 harm no, kxV [1/m], BoBrhoV [1/m], kxH, BoBrhoH, phi
                    ...

     cavity:	 cavity voltage/beam energy [eV], omega/c, beam energy [eV],
                 phi

     thin kick:	 hor., ver. displacement, roll angle (total)
		 no of nonzero multipole coeff.
		 n, b , a
		     n   n
		    ...
*/


void prt_name(FILE *fp, const int i,
	     const int type, const int method, const int n_step)
{
  fprintf(fp, "%-15s %4d %4d %4d\n",
	  elem[i].Name, elem[i].Fnum, elem[i].Knum, i);
  fprintf(fp, " %3d %3d %3d\n", type, method, n_step);
  fprintf(fp, " %23.16e %23.16e %23.16e %23.16e\n",
	  elem[i].max_ampl[X_][0], elem[i].max_ampl[X_][1],
	  elem[i].max_ampl[Y_][0], elem[i].max_ampl[Y_][1]);
}


void prt_bn(FILE *fp, const int n_design, const double bn[], const double an[],
	    const int order)
{
  int  i, nmpole;
  
  nmpole = 0;
  for (i = 0; i < order; i++)
    if ((bn[i] != 0.0) || (an[i] != 0.0)) nmpole++;
  fprintf(fp, "  %2d %2d\n", nmpole, n_design);
  for (i = 0; i < order; i++) {
    if ((bn[i] != 0.0) || (an[i] != 0.0))
      fprintf(fp, "%3d %23.16e %23.16e\n", i+1, bn[i], an[i]);
  }
}


void prt_mfile(const char file_name[])
{
  int   i, j;
  FILE  *outf;

  outf = file_write(file_name);

  for (i = 0; i < n_elem; i++) {
    switch (elem[i].kind) {
    case Drift:
      prt_name(outf, i, Drift, 0, 0);
      fprintf(outf, " %23.16e\n", elem[i].L);
      break;
    case Mpole:
      if (elem[i].L != 0.0) {
	prt_name(outf, i, Mpole, elem[i].mpole->method,
		 elem[i].mpole->n_step);
	fprintf(outf, " %23.16e %23.16e %23.16e %23.16e\n",
		elem[i].dx[X_], elem[i].dx[Y_],
		elem[i].mpole->droll_par,
		elem[i].mpole->droll_sys
		+elem[i].mpole->droll_rms*elem[i].mpole->droll_rnd);
	fprintf(outf, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		elem[i].L, elem[i].mpole->h_bend,
		elem[i].mpole->edge1, elem[i].mpole->edge2,
		elem[i].mpole->gap);
	prt_bn(outf, elem[i].mpole->n_design, elem[i].mpole->bn,
	       elem[i].mpole->an, elem[i].mpole->order);
      } else {
	prt_name(outf, i, Thinkick, elem[i].mpole->method,
		elem[i].mpole->n_step);
	fprintf(outf, " %23.16e %23.16e %23.16e\n",
		elem[i].dx[X_], elem[i].dx[Y_],
		elem[i].mpole->droll_sys
		+elem[i].mpole->droll_rms*elem[i].mpole->droll_rnd);
	prt_bn(outf, elem[i].mpole->n_design, elem[i].mpole->bn,
	       elem[i].mpole->an, elem[i].mpole->order);
      }
      break;
    case Wiggler:
      prt_name(outf, i, Wiggler, elem[i].wiggler->method,
	       elem[i].wiggler->n_step);
      fprintf(outf, " %23.16e %23.16e\n",
	      elem[i].L, elem[i].wiggler->lambda);
      fprintf(outf, "%2d\n", elem[i].wiggler->n_harm);
      for (j = 0; j < elem[i].wiggler->n_harm; j++) {
	fprintf(outf, "%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		elem[i].wiggler->harm[j],
		elem[i].wiggler->kxV[j], elem[i].wiggler->BoBrhoV[j],
		elem[i].wiggler->kxH[j], elem[i].wiggler->BoBrhoH[j],
		elem[i].wiggler->phi[j]);
      }
      break;
    case Cavity:
      prt_name(outf, i, Cavity, 0, 0);
      fprintf(outf, " %23.16e %23.16e %d %23.16e %23.16e\n",
	      elem[i].cavity->V_rf/(1e9*E0),
	      2.0*M_PI*elem[i].cavity->f_rf/clight, elem[i].cavity->h_rf,
	      1e9*E0, elem[i].cavity->phi_rf);
      break;
    case Marker:
      prt_name(outf, i, Marker, 0, 0);
      break;
    case Kick_map:
      prt_name(outf, i, Kick_map, elem[i].kick_map->method,
	       elem[i].kick_map->n_step);
      if (elem[i].kick_map->order == 1)
	fprintf(outf, " %3.1lf %1d %s\n",
		elem[i].kick_map->scl, 1, elem[i].kick_map->file_name);
      else
	fprintf(outf, " %3.1lf %1d %s\n",
		elem[i].kick_map->scl, 2, elem[i].kick_map->file_name);
      break;
    case Map_:
      prt_name(outf, i, Map_, 0, 0);
      fprintf(outf, " %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e"
	      " %23.16e %23.16e\n",
	      elem[i].map->dnu[X_],elem[i].map->dnu[Y_],
	      elem[i].map->alpha[X_],elem[i].map->alpha[Y_],
	      elem[i].map->beta[X_],elem[i].map->beta[Y_],
	      elem[i].map->eta_x,elem[i].map->etap_x);
      break;
    default:
      printf("prt_mfile: unknown type %d\n", elem[i].kind);
      exit(1);
      break;
    }
  }

  fclose(outf);
}
