#include "sxsdk.cxx"

#include "boost/scoped_array.hpp"
using namespace boost;

#define ORIGINAL_KNOT	1
#define USE_NUMERICAL_RECIPES	1
#include "numerical_recipes.hpp"

#define NUBSApproximate_id 0x30670000
#define NUBS_AUTOSMOOTH_id	uuid_class(NUBSApproximate_id)

static 	const sx::uuid_class APPROXIMATE_UUID("0A2C8F87-A7D2-4ca3-908F-7B25C2BCBB85");

namespace {
	class approximation_curve_component : public sxsdk::modifier_interface {
	public:
		approximation_curve_component(sxsdk::shade_interface *shade) : shade(shade) { }
	private:
		sxsdk::shade_interface *const shade;
		virtual int get_shade_version () const { return SHADE_BUILD_NUMBER; }
		virtual sx::uuid_class get_uuid (void *) { return APPROXIMATE_UUID; }
		virtual bool undoable (void *) const { return false; }
		virtual bool modify (sxsdk::scene_interface *scene, void *);

		int m, pmax;
		sxsdk::vec3 *vQ, *vP;
		double *t;
		int *N1;	//[7+DMAX-2];
		struct _N2 {
			int sw;
			double _t_, _1_;
		} *N2;	//[6+DMAX-2];
		struct _N3 {
			int sw;
			double _tt_, _t_, _1_;
		} *N3;	//[5+DMAX-2];
		struct _N4 {
			double _ttt_, _tt_, _t_, _1_;
		} **N4;	//[DMAX][4+DMAX-2];

		struct q_coef_round_order {
			double *q_coef;
			sxsdk::vec3 vC;
		};
		sxsdk::mat4 afN;
		sxsdk::mat4 mtx_wu;
		sxsdk::mat4 u2w;
		void Nmkr( int n ) {
			double _t_, _1_, val;
			int i;

			m = 2;
			N2[0].sw = N2[1].sw = 0;
			N2[3+pmax].sw = N2[2+pmax].sw = 0;
			for( i=2; i<2+pmax; i++ ) {
				N2[i].sw = 0;
				N2[i]._t_ = N2[i]._1_ = 0.0;
				if ( N1[i] ) {
					N2[i].sw = 1;
					val = t[i+m-1] - t[i];
					if ( 0.0 < val ) {
						N2[i]._t_ += 1.0 / val;
						N2[i]._1_ += -t[i] / val;
					}
				}
				if ( N1[i+1] ) {
					N2[i].sw = 1;
					val = t[i+m] - t[i+1];
					if ( 0.0 < val ) {
						N2[i]._1_ += t[i+m] / val;
						N2[i]._t_ += -1.0 / val;
					}
				}
			}
			m = 3;
			N3[0].sw = N3[5+pmax-3].sw = 0;
			for( i=1; i<2+pmax; i++ ) {
				N3[i].sw = 0;
				N3[i]._tt_ = N3[i]._t_ = N3[i]._1_ = 0.0;
				if ( N2[i].sw ) {
					N3[i].sw = 1;
					val = t[i+m-1] - t[i];
					if ( 0.0 < val ) {
						_t_ = 1.0 / val;
						_1_ = -t[i] / val;
						N3[i]._tt_ += _t_ * N2[i]._t_;
						N3[i]._t_ += _1_ * N2[i]._t_ + _t_ * N2[i]._1_;
						N3[i]._1_ += _1_ * N2[i]._1_;
					}
				}
				if ( N2[i+1].sw ) {
					N3[i].sw = 1;
					val = t[i+m] - t[i+1];
					if ( 0.0 < val ) {
						_1_ = t[i+m] / val;
						_t_ = -1.0 / val;
						N3[i]._tt_ += _t_ * N2[i+1]._t_;
						N3[i]._t_ += _1_ * N2[i+1]._t_ + _t_ * N2[i+1]._1_;
						N3[i]._1_ += _1_ * N2[i+1]._1_;
					}
				}
			}
			m = 4;
			for( i=0; i<2+pmax; i++ ) {
				N4[n][i]._ttt_ = N4[n][i]._tt_ = N4[n][i]._t_ = N4[n][i]._1_ = 0.0;
				if ( N3[i].sw ) {
					val = t[i+m-1] - t[i];
					if ( 0.0 < val ) {
						_t_ = 1.0 / val;
						_1_ = -t[i] / val;
						N4[n][i]._ttt_ += _t_ * N3[i]._tt_;
						N4[n][i]._tt_ += _t_ * N3[i]._t_ + _1_ * N3[i]._tt_;
						N4[n][i]._t_ += _t_ * N3[i]._1_ + _1_ * N3[i]._t_;
						N4[n][i]._1_ += _1_ * N3[i]._1_;
					}
				}
				if ( N3[i+1].sw ) {
					val = t[i+m] - t[i+1];
					if ( 0.0 < val ) {
						_1_ = t[i+m] / val;
						_t_ = -1.0 / val;
						N4[n][i]._ttt_ += _t_ * N3[i+1]._tt_;
						N4[n][i]._tt_ += _t_ * N3[i+1]._t_ + _1_ * N3[i+1]._tt_;
						N4[n][i]._t_ += _t_ * N3[i+1]._1_ + _1_ * N3[i+1]._t_;
						N4[n][i]._1_ += _1_ * N3[i+1]._1_;
					}
				}
			}
		}
		double Nj4( int i, double t, int n ) {
			return( ( (N4[n][i]._ttt_*t + N4[n][i]._tt_)*t + N4[n][i]._t_ )*t + N4[n][i]._1_ );
		}

		sxsdk::vec3 NUBS_Spline_point( double tt, int i ) {
			sxsdk::vec3 V;
			V.x = Nj4(i,   tt, i)*vQ[i].x +
					Nj4(i+1, tt, i)*vQ[i+1].x +
					Nj4(i+2, tt, i)*vQ[i+2].x +
					Nj4(i+3, tt, i)*vQ[i+3].x;
			V.y = Nj4(i,   tt, i)*vQ[i].y +
					Nj4(i+1, tt, i)*vQ[i+1].y +
					Nj4(i+2, tt, i)*vQ[i+2].y +
					Nj4(i+3, tt, i)*vQ[i+3].y;
			V.z = Nj4(i,   tt, i)*vQ[i].z +
					Nj4(i+1, tt, i)*vQ[i+1].z +
					Nj4(i+2, tt, i)*vQ[i+2].z +
					Nj4(i+3, tt, i)*vQ[i+3].z;
			return( V );
		}
		double gradientN( int i, double t, int n ) {	// ( 3.0f * N4[n][i]._ttt_ * t + 2.0f * N4[n][i]._tt_ ) * t + N4[n][i]._t_
			return( ( (N4[n][i]._ttt_+N4[n][i]._ttt_+N4[n][i]._ttt_) * t + N4[n][i]._tt_ + N4[n][i]._tt_ ) * t + N4[n][i]._t_ );
		}
		sxsdk::vec3 NUBSgradient( double tt, int i ) {
			sxsdk::vec3 V;
			V.x = gradientN(i,   t[3+i]+tt, i)*vQ[i].x +
					gradientN(i+1, t[3+i]+tt, i)*vQ[i+1].x +
					gradientN(i+2, t[3+i]+tt, i)*vQ[i+2].x +
					gradientN(i+3, t[3+i]+tt, i)*vQ[i+3].x;
			V.y = gradientN(i,   t[3+i]+tt, i)*vQ[i].y +
					gradientN(i+1, t[3+i]+tt, i)*vQ[i+1].y +
					gradientN(i+2, t[3+i]+tt, i)*vQ[i+2].y +
					gradientN(i+3, t[3+i]+tt, i)*vQ[i+3].y;
			V.z = gradientN(i,   t[3+i]+tt, i)*vQ[i].z +
					gradientN(i+1, t[3+i]+tt, i)*vQ[i+1].z +
					gradientN(i+2, t[3+i]+tt, i)*vQ[i+2].z +
					gradientN(i+3, t[3+i]+tt, i)*vQ[i+3].z;
			return( V );
		}
		sxsdk::vec3 NUBSgradient_( double tt, int i ) {
			sxsdk::vec3 V;
			V.x = gradientN(i,   t[4+i]+tt, i)*vQ[i].x +
					gradientN(i+1, t[4+i]+tt, i)*vQ[i+1].x +
					gradientN(i+2, t[4+i]+tt, i)*vQ[i+2].x +
					gradientN(i+3, t[4+i]+tt, i)*vQ[i+3].x;
			V.y = gradientN(i,   t[4+i]+tt, i)*vQ[i].y +
					gradientN(i+1, t[4+i]+tt, i)*vQ[i+1].y +
					gradientN(i+2, t[4+i]+tt, i)*vQ[i+2].y +
					gradientN(i+3, t[4+i]+tt, i)*vQ[i+3].y;
			V.z = gradientN(i,   t[4+i]+tt, i)*vQ[i].z +
					gradientN(i+1, t[4+i]+tt, i)*vQ[i+1].z +
					gradientN(i+2, t[4+i]+tt, i)*vQ[i+2].z +
					gradientN(i+3, t[4+i]+tt, i)*vQ[i+3].z;
			return( V );
		}

		void convert_NUBS_to_Bezier( sxsdk::line_class &l_imp, int ls=-1, int le=-1, int line_pmax=0, int dir=0 ) {
			int i, li, ib;
			sxsdk::vec3 vGl, vGl2;
			sxsdk::vec3 bef_ihnd, bef_ohnd, bef2_ihnd, bef2_ohnd;
			sxsdk::vec3 v3d;
			sxsdk::vec3 v3;
			bool has_bef_ihnd=false;
			if ( ls == -1 )	ls = 0;
			if ( le == -1 )	le = pmax-1;
			int ps = 0;
			//if ( dir )	ps = gContinuous_points_i;
			int ls_bak = 0, le_bak;
			if ( le < ls ) {
				ps = line_pmax - ls;
				ls_bak = ls;
				ls = 0;
			}
			le++;
			int le_1 = le - 1;
			int le_2 = le - 2;

			sxsdk::mat4 w2l = l_imp.get_local_to_world_matrix().inv();
			sxsdk::mat4 uwl = u2w * w2l;
			uwl = sxsdk::mat4::identity;

			for( i=ib=ps, li=ls; li<le_1; i++,li++,ib++ ) {
				if ( i<pmax-1 ) {
					vGl = NUBSgradient( 0.0, i );
					vGl.x *= ( t[4+i] - t[3+i] );
					vGl.y *= ( t[4+i] - t[3+i] );
					vGl.z *= ( t[4+i] - t[3+i] );
					vGl.x /= 3.0;
					vGl.y /= 3.0;
					vGl.z /= 3.0;
					vGl += vP[i];

					if ( li < le_2 )
						vGl2 = NUBSgradient( 0.0, i+1 );
					else
						vGl2 = NUBSgradient_( 0.0, i );
					vGl2 = -vGl2;
					vGl2.x *= ( t[4+i] - t[3+i] );
					vGl2.y *= ( t[4+i] - t[3+i] );
					vGl2.z *= ( t[4+i] - t[3+i] );
					vGl2.x /= 3.0;
					vGl2.y /= 3.0;
					vGl2.z /= 3.0;
					vGl2 += vP[i+1];
					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					#ifdef AUTO_SMOOTH_DBG
					{logging_class dbg("shade_dbg.fil"); if (dbg) fprintf(dbg,"vGl:%f %f %f, vGl2:%f %f %f\n", vGl.p[0], vGl.p[1], vGl.p[2], vGl2.p[0], vGl2.p[1], vGl2.p[2]); }
					#endif
					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					if ( has_bef_ihnd ) {
						bef2_ohnd = vGl;
						bef2_ihnd = vGl2;
						l_imp.set_outhandle_deprecated( li, bef2_ohnd * uwl );	l_imp.set_handle_linked_deprecated( li, true );
						l_imp.set_inhandle_deprecated( li, bef_ihnd * uwl );	l_imp.set_handle_linked_deprecated( li, true );
					} else {
						l_imp.set_outhandle_deprecated( li, vGl * uwl );	l_imp.set_handle_linked_deprecated( li, true );
					}
					bef_ohnd = vGl;
					bef_ihnd = vGl2;
					has_bef_ihnd = true;
				} else {
					v3d = vP[ib] ;//* uws;
					v3 = v3d;
					l_imp.set_outhandle_deprecated( li, v3 * uwl );	l_imp.set_handle_linked_deprecated( li, true );
					v3d = vP[ib+1] ;//* uws;
					v3 = v3d;
					l_imp.set_inhandle_deprecated( li+1, v3 * uwl );	l_imp.set_handle_linked_deprecated( li, true );
				}
			}
			if ( i-1<pmax-1 ) {
				if ( li<le ) {
					l_imp.set_inhandle_deprecated( li, bef_ihnd * uwl );	l_imp.set_handle_linked_deprecated( li, true );
				}
			}
		}
	};

}
bool approximation_curve_component::modify (sxsdk::scene_interface *scene, void *) {
	try {
		int i=0, j;

		int num = scene->get_active_shapes(0);
		scoped_array<sxsdk::shape_class*> sel_shapes(new  sxsdk::shape_class *[ num ]);
		scene->get_active_shapes( sel_shapes.get() );
		
		// \en Get a stream for accessing custom attributes of the shape.\enden \ja 形状のカスタム属性にアクセスするためのストリームを取得する。 \endja 
		compointer<sxsdk::stream_interface> stream(sel_shapes[i]->get_attribute_stream_interface_with_uuid(APPROXIMATE_UUID));

		if ( sel_shapes[i]->get_type() == sxsdk::enums::line ) {
			sxsdk::line_class &line = sel_shapes[i]->get_line();
			int line_pmax = line.get_total_number_of_control_points();
			int appried_num = 0;
			std::vector<int> appried_at;
			std::vector<float> appried_weights;

			bool jp = scene->is_japanese_mode();
			compointer<sxsdk::strings_interface> strings(scene->create_strings_interface("text"));

			// \en Create a dialog box inquirering a numerical value.\enden \ja 数値を1個だけ入力するダイアログボックスを作成する。 \endja 
			compointer<sxsdk::dialog_interface> dialog(shade->create_dialog_interface());
			dialog->set_title( strings->gettext("ApproximationCurve_title") );
			char buf[64];
			sprintf(buf, "%s (1~%d)", strings->gettext("ApproximationCurve_seg_num"), line_pmax-3);
			dialog->append_item(sxsdk::dialog_interface::int_type, buf);
			//dialog->append_item(sxsdk::dialog_interface::int_type, "appried weighted points number");
			dialog->append_item( sxsdk::dialog_interface::float_type,  strings->gettext("ApproximationCurve_weight")  );

			if (stream) {
				// \en If the custom attributes exists, read a floating point value from the attribute stream and set it to the dialog box.\enden \ja カスタム属性が存在する場合、そこから浮動小数点数１個を読み込み、ダイアログボックスにその値を設定する。 \endja 
				float ver;
				int seg;
				stream->read_float(ver);
				stream->read_int(seg);
				stream->read_int(appried_num);
				int cur_appried_at;
				float cur_appried_weight;
				for( int i=0; i<appried_num; i++ ) {
					stream->read_int(cur_appried_at);
					stream->read_float(cur_appried_weight);
					appried_at.push_back(cur_appried_at);
					appried_weights.push_back(cur_appried_weight);
				}

				dialog->set_int_property_value(0, seg);
				//dialog->set_int_property_value(1, appried_num);

				int cur_appried_num = line.get_active_control_points(0);
				if ( cur_appried_num == 1 ) {
					line.get_active_control_points(&cur_appried_at);
					bool mutched = false;
					for( int j=0; j<appried_num; j++ ) {
						if ( 0 <= cur_appried_at && cur_appried_at < line_pmax ) {
							if ( cur_appried_at == appried_at[j] ) {
								cur_appried_weight = appried_weights[j];
								mutched = true;
								break;
							}
						}
					}
					if ( mutched )	dialog->set_float_property_value(1, cur_appried_weight);
					else			dialog->set_float_property_value(1, 1.0);
				} else			dialog->set_float_property_value(1, 1.0);
			} else {
				dialog->set_int_property_value(0, line_pmax/2);
				//dialog->set_int_property_value(1, 0);
				dialog->set_float_property_value(1, 1.0);
			}
			if (dialog->ask()) { // \en Display the dialog box.\enden \ja ダイアログボックスを表示する。 \endja 
				int seg, dummy;
				float current_weight;
				dialog->get_int_property_value(0, seg);
				//dialog->get_int_property_value(1, dummy);
				dialog->get_float_property_value(1, current_weight);
				
				compointer<sxsdk::stream_interface> stream(sel_shapes[i]->create_attribute_stream_interface_with_uuid(APPROXIMATE_UUID));
				stream->write_float(1.0f);	// version num
				stream->write_int(seg);
				stream->set_label("approximation_curve");

				int new_appried_num = line.get_active_control_points(0);
				if ( 0 < new_appried_num ) {
					int *new_appried_at = (int *)malloc(sizeof(int)*new_appried_num);
					line.get_active_control_points(new_appried_at);

					if ( 0 < appried_num ) {
						for( int i=0; i<new_appried_num; i++ ) {
							bool mutched = false;
							for( int j=0; j<appried_num; j++ ) {
								if ( 0 <= new_appried_at[i] && new_appried_at[i] < line_pmax ) {
									if ( new_appried_at[i] == appried_at[j] ) {
										appried_weights[j] = current_weight;
										mutched = true;
										break;
									}
								}
							}
							if ( mutched == false ) {
								appried_at.push_back(new_appried_at[i]);
								appried_weights.push_back(current_weight);
								appried_num++;
							}
						}
					} else {
						for( int i=0; i<new_appried_num; i++ ) {
							if ( 0 <= new_appried_at[i] && new_appried_at[i] < line_pmax ) {
								appried_at.push_back(new_appried_at[i]);
								appried_weights.push_back(current_weight);
							}
						}
						appried_num = new_appried_num;
					}
				}
				stream->write_int(appried_num);
				for( int i=0; i<appried_num; i++ ) {
					stream->write_int(appried_at[i]);
					stream->write_float(appried_weights[i]);
				}

				sxsdk::mat4 l2w = line.get_local_to_world_matrix();
				sxsdk::mat4 w2u = scene->get_world_to_user_matrix();
				sxsdk::mat4 lwu = l2w * w2u;
				u2w = scene->get_user_to_world_matrix();

				int nSamples = line.get_total_number_of_control_points();
				int nCtrlPnts;
				stream->set_pointer(0);
				float ver;
				stream->read_float(ver);
				stream->read_int(nCtrlPnts);
				nCtrlPnts += 3;
				if ( nCtrlPnts < 4 )	nCtrlPnts = 4;
				if ( nSamples < nCtrlPnts ) nCtrlPnts = nSamples;

				sxsdk::vec3 vTmp;

				vP = (sxsdk::vec3 *)malloc( sizeof(sxsdk::vec3) * nSamples );
				for( i=0; i<nSamples; i++ ) {
					vP[i] = line.get_anchor_point_deprecated(i) * lwu;
				}
				int rank = 4;
				double *_t = (double *)malloc( sizeof(double) * nSamples );
				double *ls = (double *)malloc( sizeof(double) * nSamples );
				ls[0] = 0.0;
				double len_sigma = 0.0;
				sxsdk::vec3 dif;
				for( i=1; i<nSamples; i++ ) {
					dif = vP[i] - vP[i-1];
					len_sigma += ls[i] = sqrt( dif.x*dif.x + dif.y*dif.y + dif.z*dif.z );
				}
				double len = 0.0;
				for( i=0; i<nSamples; i++ ) {
					len += ls[i];
					_t[i] = len / len_sigma;
				}

				t = (double *)malloc( sizeof(double) * (nCtrlPnts+rank+3) );
				#if ORIGINAL_KNOT==0
				for( i=0; i<rank; i++ )	t[i] = 0.0;
				int j=1;
				double Alpha, Delta = (double)nSamples / (double)(nCtrlPnts - rank + 1);
				double multiplier = 1.0 / (double)(nCtrlPnts - 3);
				for( i=j+rank-1; i<nCtrlPnts; j++,i++ ) {
					#if 1
					int ii = (int)( (double)j * Delta );
					Alpha = (double)j * Delta - (double)ii;
					t[i] = (1.0 - Alpha) * _t[ii-1] + Alpha * _t[ii];
					#else
					t[i] = (double)j * multiplier;
					#endif
				}
				for( ; i<nCtrlPnts+rank+3; i++ )	t[i] = 1.0;
				#else
				for( i=0; i<rank; i++ )	t[i] = 0.0;

				int segment_num = nCtrlPnts - 3;
				int undef_knot_num = segment_num - 1;
				int influ_num = Round((double)nSamples * ( 2.0 / (double)(nCtrlPnts-2) ), 0);
				int step_num = nSamples / segment_num;
				sxsdk::vec3 *knotP = (sxsdk::vec3 *)malloc( sizeof(sxsdk::vec3) * (nCtrlPnts-2) );
				knotP[0] = vP[0];
				knotP[segment_num] = vP[nSamples-1];

				int influ_half = influ_num / 2;
				int shift = step_num;
				if ( shift < influ_half )	shift = influ_half;
				int idx;
				for( idx=0; idx<undef_knot_num/2; idx++, shift+=step_num ) {
					knotP[idx+1] = sxsdk::vec3(0.0, 0.0, 0.0);
					for( i=0; i<influ_num; i++ ) {
						knotP[idx+1] += vP[shift-influ_half+i];
					}
					knotP[idx+1].x /= (double)influ_num;
					knotP[idx+1].y /= (double)influ_num;
					knotP[idx+1].z /= (double)influ_num;
				}
				shift = nSamples-1 - step_num;
				if ( nSamples-1 < shift+influ_half )	shift = nSamples-1 - influ_half;
				for( idx=undef_knot_num-1; undef_knot_num/2<=idx; idx--, shift-=step_num ) {
					knotP[idx+1] = sxsdk::vec3(0.0, 0.0, 0.0);
					for( i=0; i<influ_num; i++ ) {
						knotP[idx+1] += vP[shift+influ_half-i];
					}
					knotP[idx+1].x /= (double)influ_num;
					knotP[idx+1].y /= (double)influ_num;
					knotP[idx+1].z /= (double)influ_num;
				}
				if ( undef_knot_num % 2 ) {
					int mid_i = nSamples / 2;
					idx = undef_knot_num / 2;
					knotP[idx+1] = sxsdk::vec3(0.0, 0.0, 0.0);
					if ( nSamples % 2 ) {
						for( i=0; i<3; i++ ) {
							knotP[idx+1] += vP[mid_i-1+i];
						}
						knotP[idx+1].x /= (double)3;
						knotP[idx+1].y /= (double)3;
						knotP[idx+1].z /= (double)3;
					} else {
						for( i=0; i<4; i++ ) {
							knotP[idx+1] += vP[mid_i-2+i];
						}
						knotP[idx+1].x /= (double)4;
						knotP[idx+1].y /= (double)4;
						knotP[idx+1].z /= (double)4;
					}
				}

				if ( 2 <= undef_knot_num ) {
					sxsdk::vec3 knotP1, knotPn_1;
					knotP1 = (knotP[0]*sxsdk::vec3(1.0) + knotP[2]*sxsdk::vec3(4.0)) / sxsdk::vec3(5.0);
					knotPn_1 = (knotP[segment_num]*sxsdk::vec3(1.0) + knotP[segment_num-2]*sxsdk::vec3(4.0)) / sxsdk::vec3(5.0);
					knotP[1] = knotP1;
					knotP[segment_num-1] = knotPn_1;
				}

				len_sigma = 0.0;
				for( i=1; i<undef_knot_num+2; i++ ) {
					dif = knotP[i] - knotP[i-1];
					len_sigma += ls[i-1] = sqrt( dif.x*dif.x + dif.y*dif.y + dif.z*dif.z );
				}
				len = 0.0;
				for( i=0; i<undef_knot_num; i++ ) {
					len += ls[i];
					t[rank+i] = len / len_sigma;
				}

				for( i=nCtrlPnts; i<nCtrlPnts+rank+3; i++ )	t[i] = 1.0;
				#endif

				free( ls );

				pmax = nCtrlPnts + 1;
				N1 = (int *)malloc( sizeof(int) * (7+nCtrlPnts-1) );	//[7+DMAX-2];
				N2 = (_N2 *)malloc( sizeof(_N2) * (6+nCtrlPnts-1) );	//[6+DMAX-2];
				N3 = (_N3 *)malloc( sizeof(_N3) * (5+nCtrlPnts-1) );	//[5+DMAX-2];
				N4 = (_N4 **)malloc( sizeof(_N4 *) * nCtrlPnts );
				for( i=0; i<nCtrlPnts; i++ ) {
					N4[i] = (_N4 *)malloc( sizeof(_N4) * (4+nCtrlPnts-1) );
					for( j=0; j<nCtrlPnts+3; j++ ) {
						N4[i][j]._1_ = 0.0;
						N4[i][j]._t_ = 0.0;
						N4[i][j]._tt_ = 0.0;
						N4[i][j]._ttt_ = 0.0;
					}
				}
				for( i=0; i<=4+nCtrlPnts; i++ ) N1[i] = 0;
				for( i=3; i<=2+nCtrlPnts; i++ )  {
					N1[i] = 1;
					Nmkr(i-3);
					N1[i] = 0;
				}

				q_coef_round_order *round_order = (q_coef_round_order *)malloc( sizeof(q_coef_round_order) * (nCtrlPnts) );
				for( i=0; i<nCtrlPnts; i++ ) {
					round_order[i].q_coef = (double *)malloc( sizeof(double) * (nCtrlPnts) );
					for( j=0; j<nCtrlPnts; j++ ) {
						round_order[i].q_coef[j] = 0.0;
					}
					round_order[i].vC = sxsdk::vec3(0.0, 0.0, 0.0);
				}

				double coef_q0, coef_q1, coef_q2, coef_q3, coef_q4, coef_q5;
				int ii, cur_seg_i;
				for( i=0; i<nSamples; i++ ) {
					float cur_weight = 1.0;
					for( int j=0; j<appried_num; j++ ) {
						if ( i == appried_at[j] ) {
							cur_weight = appried_weights[j];
							break;
						}
					}
					// get sample point's curve_segment
					for( ii=rank; ii<nCtrlPnts; ii++ ) {
						if ( t[ii-1] < _t[i] && _t[i] <= t[ii] )
							break;
					}
					cur_seg_i = ii - rank;	//
					if ( i == 0 )	cur_seg_i = 0;

					coef_q1 = Nj4( cur_seg_i,   _t[i], cur_seg_i );
					coef_q2 = Nj4( cur_seg_i+1, _t[i], cur_seg_i );
					coef_q3 = Nj4( cur_seg_i+2, _t[i], cur_seg_i );
					coef_q4 = Nj4( cur_seg_i+3, _t[i], cur_seg_i );

					// make round 4 kind
					if ( i < nSamples-1 ) {
						round_order[cur_seg_i  ].q_coef[cur_seg_i  ] += coef_q1 * coef_q1 * 2.0 * cur_weight;	//
						round_order[cur_seg_i  ].q_coef[cur_seg_i+1] += coef_q1 * coef_q2 * 2.0 * cur_weight;
						round_order[cur_seg_i  ].q_coef[cur_seg_i+2] += coef_q1 * coef_q3 * 2.0 * cur_weight;
						round_order[cur_seg_i  ].q_coef[cur_seg_i+3] += coef_q1 * coef_q4 * 2.0 * cur_weight;
						round_order[cur_seg_i  ].vC.x += vP[i].x * -(coef_q1 + coef_q1) * cur_weight;
						round_order[cur_seg_i  ].vC.y += vP[i].y * -(coef_q1 + coef_q1) * cur_weight;
						round_order[cur_seg_i  ].vC.z += vP[i].z * -(coef_q1 + coef_q1) * cur_weight;
					}

					round_order[cur_seg_i+1].q_coef[cur_seg_i  ] += coef_q2 * coef_q1 * 2.0 * cur_weight;
					round_order[cur_seg_i+1].q_coef[cur_seg_i+1] += coef_q2 * coef_q2 * 2.0 * cur_weight;	//
					round_order[cur_seg_i+1].q_coef[cur_seg_i+2] += coef_q2 * coef_q3 * 2.0 * cur_weight;
					round_order[cur_seg_i+1].q_coef[cur_seg_i+3] += coef_q2 * coef_q4 * 2.0 * cur_weight;
					round_order[cur_seg_i+1].vC.x += vP[i].x * -(coef_q2 + coef_q2) * cur_weight;
					round_order[cur_seg_i+1].vC.y += vP[i].y * -(coef_q2 + coef_q2) * cur_weight;
					round_order[cur_seg_i+1].vC.z += vP[i].z * -(coef_q2 + coef_q2) * cur_weight;

					round_order[cur_seg_i+2].q_coef[cur_seg_i  ] += coef_q3 * coef_q1 * 2.0 * cur_weight;
					round_order[cur_seg_i+2].q_coef[cur_seg_i+1] += coef_q3 * coef_q2 * 2.0 * cur_weight;
					round_order[cur_seg_i+2].q_coef[cur_seg_i+2] += coef_q3 * coef_q3 * 2.0 * cur_weight;	//
					round_order[cur_seg_i+2].q_coef[cur_seg_i+3] += coef_q3 * coef_q4 * 2.0 * cur_weight;
					round_order[cur_seg_i+2].vC.x += vP[i].x * -(coef_q3 + coef_q3) * cur_weight;
					round_order[cur_seg_i+2].vC.y += vP[i].y * -(coef_q3 + coef_q3) * cur_weight;
					round_order[cur_seg_i+2].vC.z += vP[i].z * -(coef_q3 + coef_q3) * cur_weight;

					if ( 0 < i ) {
						round_order[cur_seg_i+3].q_coef[cur_seg_i  ] += coef_q4 * coef_q1 * 2.0 * cur_weight;
						round_order[cur_seg_i+3].q_coef[cur_seg_i+1] += coef_q4 * coef_q2 * 2.0 * cur_weight;
						round_order[cur_seg_i+3].q_coef[cur_seg_i+2] += coef_q4 * coef_q3 * 2.0 * cur_weight;
						round_order[cur_seg_i+3].q_coef[cur_seg_i+3] += coef_q4 * coef_q4 * 2.0 * cur_weight;	//
						round_order[cur_seg_i+3].vC.x += vP[i].x * -(coef_q4 + coef_q4) * cur_weight;
						round_order[cur_seg_i+3].vC.y += vP[i].y * -(coef_q4 + coef_q4) * cur_weight;
						round_order[cur_seg_i+3].vC.z += vP[i].z * -(coef_q4 + coef_q4) * cur_weight;
					}
				}

				vQ = new sxsdk::vec3[nCtrlPnts];

				#if USE_NUMERICAL_RECIPES
				NUMERICAL_RECIPES nr;
				matrix4x4 Left_mtx, invLeft_mtx;
				Left_mtx = new_matrix4x4( nCtrlPnts, nCtrlPnts );
				invLeft_mtx = new_matrix4x4( nCtrlPnts, nCtrlPnts );

				for( i=0; i<nCtrlPnts; i++ ) {
					for( j=0; j<nCtrlPnts; j++ ) {
						Left_mtx[i][j] = round_order[i].q_coef[j];
					}
				}

				nr.invmtx( nCtrlPnts, Left_mtx, invLeft_mtx );

				#if 0
				matrix4x4 X_mtx;
				X_mtx = new_matrix4x4( nCtrlPnts, nCtrlPnts );
				for( i=0; i<nCtrlPnts; i++ ) {
					for( j=0; j<nCtrlPnts; j++ ) {
						X_mtx[i][j] = 0.0;
						for( int k=0; k<nCtrlPnts; k++ )	X_mtx[i][j] += Left_mtx[i][k] * invLeft_mtx[k][j];
					}
				}
				free_matrix4x4(X_mtx);
				#endif

				for( i=0; i<nCtrlPnts; i++ ) {
					vQ[i].x = 0.0;
					vQ[i].y = 0.0;
					vQ[i].z = 0.0;
					for( j=0; j<nCtrlPnts; j++ ) {
						vQ[i].x += invLeft_mtx[i][j] * -round_order[j].vC.x;
						vQ[i].y += invLeft_mtx[i][j] * -round_order[j].vC.y;
						vQ[i].z += invLeft_mtx[i][j] * -round_order[j].vC.z;
					}
				}

				free_matrix4x4(Left_mtx);
				free_matrix4x4(invLeft_mtx);
				#endif

				scene->begin_creating();

				#if 0
				scene->begin_line("knot polygon",false);
				for ( i = 0; i < nCtrlPnts-2; ++i)
				{
					scene->append_point( knotP[i] );
				}
				scene->end_line();

				scene->begin_line("Control Polygon",false);
				for ( i = 0; i < nCtrlPnts; ++i)
				{
					scene->append_point( vQ[i] );
				}
				scene->end_line();
				#endif

				free( vP );	vP = 0;
				vP = (sxsdk::vec3 *)malloc( sizeof(sxsdk::vec3) * nCtrlPnts );
				scene->begin_line(strings->gettext("label"),false);
				scene->append_point( vQ[0] );
				int ent_i = 0;
				vP[ent_i++] = vQ[0];
				for( i=rank, j=0; i<nCtrlPnts; i++, j++ ) {
					vTmp = NUBS_Spline_point(t[i], j);
					scene->append_point( vTmp );
					vP[ent_i++] = vTmp;
				}
				scene->append_point( vQ[nCtrlPnts-1] );
				vP[ent_i++] = vQ[nCtrlPnts-1];
				scene->end_line();
				pmax = ent_i;
				sxsdk::shape_class &cur_shape = scene->first_active_shape();
				sxsdk::line_class &guid_line = cur_shape.get_line();
				convert_NUBS_to_Bezier(guid_line);

				#if 0
				scene->begin_line("approximated polygon",false);
				scene->append_point( vQ[0] );
				for( i=1; i<nSamples-1; i++ ) {
					// get sample point's curve_segment
					for( ii=rank; ii<rank+nCtrlPnts-3; ii++ ) {
						if ( t[ii-1] < _t[i] && _t[i] <= t[ii] )
							break;
					}
					cur_seg_i = ii - rank;	//

					vTmp = NUBS_Spline_point(_t[i], cur_seg_i);
					scene->append_point( vTmp );
				}
				scene->append_point( vQ[nCtrlPnts-1] );
				scene->end_line();
				#endif

				sxsdk::mat4 afI = sxsdk::mat4::identity;
				afI = line.get_transformation();
				compointer<sxsdk::shape_interface> s( scene->get_insertion_point_shape_interface() );
				s->transform( lwu.inv() * afI );

				scene->end_creating();

				free( _t );
				free( t );
				#if ORIGINAL_KNOT
				free( knotP );
				#endif
				free( N1 );
				free( N2 );
				free( N3 );
				for( i=0; i<nCtrlPnts; i++ ) {
					free( N4[i] );
				}
				free( N4 );
				for( i=0; i<nCtrlPnts; i++ ) {
					free( round_order[i].q_coef );
				}
				free( round_order );
				delete[] vQ;
				free( vP );
			}
		}
	} catch (...) { return false; }
	return true;
}

SXPLUGINNAMESPACEBEGIN(test_modifier)
	SXPLUGINENTRY void __stdcall create_interface (const IID &iid, int i, void **p, sxsdk::shade_interface *shade, void *) {
		approximation_curve_component *u = new approximation_curve_component(shade);
		u->AddRef();
		*p = (void *)u;
	};
	SXPLUGINENTRY int __stdcall has_interface (const IID &iid, sxsdk::shade_interface *shade) {
		if (iid == modifier_iid) return 1;
		return 0;
	}
	SXPLUGINENTRY const char * __stdcall get_name (const IID &iid, int i, sxsdk::shade_interface *shade, void *) {
		compointer<sxsdk::strings_interface> strings(shade->create_strings_interface("text"));
		return strings->gettext("ApproximationCurve_title");
	}
	SXPLUGINENTRY sx::uuid_class __stdcall get_uuid (const IID &iid, int i, void *) {
		return APPROXIMATE_UUID;
	}
	SXPLUGINENTRY const char * __stdcall get_group_name (const IID &iid, sxsdk::shade_interface *shade, void *) {
		return 0;
	}
	SXPLUGINENTRY bool __stdcall is_resident (const IID &iid, int i, void *) {
		return false;
	}
SXPLUGINNAMESPACEEND
