# This file was automatically created by FeynRules 2.3.19
# Mathematica version: 10.1.0  for Mac OS X x86 (64-bit) (March 24, 2015)
# Date: Wed 30 Mar 2016 12:59:16


from object_library import all_vertices, all_CTvertices, Vertex, CTVertex
import particles as P
import CT_couplings as C
import lorentz as L


V_1 = CTVertex(name = 'V_1',
               type = 'R2',
               particles = [ P.g, P.g, P.g ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVV1, L.VVV2, L.VVV4, L.VVV6, L.VVV7, L.VVV8 ],
               loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ] ],
               couplings = {(0,0,0):C.R2GC_722_258,(0,0,1):C.R2GC_722_259,(0,1,0):C.R2GC_725_260,(0,1,1):C.R2GC_725_261,(0,2,0):C.R2GC_725_260,(0,2,1):C.R2GC_725_261,(0,3,0):C.R2GC_722_258,(0,3,1):C.R2GC_722_259,(0,4,0):C.R2GC_722_258,(0,4,1):C.R2GC_722_259,(0,5,0):C.R2GC_725_260,(0,5,1):C.R2GC_725_261})

V_2 = CTVertex(name = 'V_2',
               type = 'R2',
               particles = [ P.g, P.g, P.g, P.g ],
               color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
               lorentz = [ L.VVVV2, L.VVVV3, L.VVVV4 ],
               loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ] ],
               couplings = {(2,0,0):C.R2GC_580_132,(2,0,1):C.R2GC_580_133,(0,0,0):C.R2GC_580_132,(0,0,1):C.R2GC_580_133,(4,0,0):C.R2GC_578_128,(4,0,1):C.R2GC_578_129,(3,0,0):C.R2GC_578_128,(3,0,1):C.R2GC_578_129,(8,0,0):C.R2GC_579_130,(8,0,1):C.R2GC_579_131,(7,0,0):C.R2GC_584_139,(7,0,1):C.R2GC_731_266,(6,0,0):C.R2GC_583_137,(6,0,1):C.R2GC_732_267,(5,0,0):C.R2GC_578_128,(5,0,1):C.R2GC_578_129,(1,0,0):C.R2GC_578_128,(1,0,1):C.R2GC_578_129,(11,0,0):C.R2GC_582_135,(11,0,1):C.R2GC_582_136,(10,0,0):C.R2GC_582_135,(10,0,1):C.R2GC_582_136,(9,0,1):C.R2GC_581_134,(2,1,0):C.R2GC_580_132,(2,1,1):C.R2GC_580_133,(0,1,0):C.R2GC_580_132,(0,1,1):C.R2GC_580_133,(4,1,0):C.R2GC_578_128,(4,1,1):C.R2GC_578_129,(3,1,0):C.R2GC_578_128,(3,1,1):C.R2GC_578_129,(8,1,0):C.R2GC_579_130,(8,1,1):C.R2GC_733_268,(6,1,0):C.R2GC_728_262,(6,1,1):C.R2GC_728_263,(7,1,0):C.R2GC_584_139,(7,1,1):C.R2GC_584_140,(5,1,0):C.R2GC_578_128,(5,1,1):C.R2GC_578_129,(1,1,0):C.R2GC_578_128,(1,1,1):C.R2GC_578_129,(11,1,0):C.R2GC_582_135,(11,1,1):C.R2GC_582_136,(10,1,0):C.R2GC_582_135,(10,1,1):C.R2GC_582_136,(9,1,1):C.R2GC_581_134,(2,2,0):C.R2GC_580_132,(2,2,1):C.R2GC_580_133,(0,2,0):C.R2GC_580_132,(0,2,1):C.R2GC_580_133,(4,2,0):C.R2GC_578_128,(4,2,1):C.R2GC_578_129,(3,2,0):C.R2GC_578_128,(3,2,1):C.R2GC_578_129,(8,2,0):C.R2GC_579_130,(8,2,1):C.R2GC_730_265,(6,2,0):C.R2GC_583_137,(6,2,1):C.R2GC_583_138,(7,2,0):C.R2GC_729_264,(7,2,1):C.R2GC_580_133,(5,2,0):C.R2GC_578_128,(5,2,1):C.R2GC_578_129,(1,2,0):C.R2GC_578_128,(1,2,1):C.R2GC_578_129,(11,2,0):C.R2GC_582_135,(11,2,1):C.R2GC_582_136,(10,2,0):C.R2GC_582_135,(10,2,1):C.R2GC_582_136,(9,2,1):C.R2GC_581_134})

V_3 = CTVertex(name = 'V_3',
               type = 'R2',
               particles = [ P.d__tilde__, P.u, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.d, P.g, P.u] ] ],
               couplings = {(0,0,0):C.R2GC_664_208,(0,1,0):C.R2GC_659_203})

V_4 = CTVertex(name = 'V_4',
               type = 'R2',
               particles = [ P.s__tilde__, P.u, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.g, P.s, P.u] ] ],
               couplings = {(0,0,0):C.R2GC_670_214,(0,1,0):C.R2GC_665_209})

V_5 = CTVertex(name = 'V_5',
               type = 'R2',
               particles = [ P.b__tilde__, P.u, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.b, P.g, P.u] ] ],
               couplings = {(0,0,0):C.R2GC_718_253,(0,1,0):C.R2GC_674_215})

V_6 = CTVertex(name = 'V_6',
               type = 'R2',
               particles = [ P.d__tilde__, P.c, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.c, P.d, P.g] ] ],
               couplings = {(0,0,0):C.R2GC_628_172,(0,1,0):C.R2GC_623_167})

V_7 = CTVertex(name = 'V_7',
               type = 'R2',
               particles = [ P.s__tilde__, P.c, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.c, P.g, P.s] ] ],
               couplings = {(0,0,0):C.R2GC_634_178,(0,1,0):C.R2GC_629_173})

V_8 = CTVertex(name = 'V_8',
               type = 'R2',
               particles = [ P.b__tilde__, P.c, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.b, P.c, P.g] ] ],
               couplings = {(0,0,0):C.R2GC_719_254,(0,1,0):C.R2GC_675_216})

V_9 = CTVertex(name = 'V_9',
               type = 'R2',
               particles = [ P.d__tilde__, P.t, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS4 ],
               loop_particles = [ [ [P.d, P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_778_304,(0,1,0):C.R2GC_737_269})

V_10 = CTVertex(name = 'V_10',
                type = 'R2',
                particles = [ P.s__tilde__, P.t, P.h__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.s, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_779_305,(0,1,0):C.R2GC_738_270})

V_11 = CTVertex(name = 'V_11',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.h__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_789_315,(0,1,0):C.R2GC_739_271})

V_12 = CTVertex(name = 'V_12',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_663_207,(0,1,0):C.R2GC_662_206})

V_13 = CTVertex(name = 'V_13',
                type = 'R2',
                particles = [ P.c__tilde__, P.d, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_627_171,(0,1,0):C.R2GC_625_169})

V_14 = CTVertex(name = 'V_14',
                type = 'R2',
                particles = [ P.t__tilde__, P.d, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_772_298,(0,1,0):C.R2GC_766_292})

V_15 = CTVertex(name = 'V_15',
                type = 'R2',
                particles = [ P.u__tilde__, P.s, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.s, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_669_213,(0,1,0):C.R2GC_668_212})

V_16 = CTVertex(name = 'V_16',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_633_177,(0,1,0):C.R2GC_631_175})

V_17 = CTVertex(name = 'V_17',
                type = 'R2',
                particles = [ P.t__tilde__, P.s, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.s, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_773_299,(0,1,0):C.R2GC_767_293})

V_18 = CTVertex(name = 'V_18',
                type = 'R2',
                particles = [ P.u__tilde__, P.b, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_707_242,(0,1,0):C.R2GC_699_234})

V_19 = CTVertex(name = 'V_19',
                type = 'R2',
                particles = [ P.c__tilde__, P.b, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_708_243,(0,1,0):C.R2GC_700_235})

V_20 = CTVertex(name = 'V_20',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_774_300,(0,1,0):C.R2GC_768_294})

V_21 = CTVertex(name = 'V_21',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_770_296,(0,1,0):C.R2GC_771_297})

V_22 = CTVertex(name = 'V_22',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_692_230})

V_23 = CTVertex(name = 'V_23',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_759_288,(0,1,0):C.R2GC_757_286})

V_24 = CTVertex(name = 'V_24',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_758_287})

V_25 = CTVertex(name = 'V_25',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_620_164,(0,1,0):C.R2GC_617_161})

V_26 = CTVertex(name = 'V_26',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_621_165,(0,1,0):C.R2GC_618_162})

V_27 = CTVertex(name = 'V_27',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_622_166,(0,1,0):C.R2GC_619_163})

V_28 = CTVertex(name = 'V_28',
                type = 'R2',
                particles = [ P.c__tilde__, P.u, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_656_200,(0,1,0):C.R2GC_647_191})

V_29 = CTVertex(name = 'V_29',
                type = 'R2',
                particles = [ P.c__tilde__, P.u, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_657_201,(0,1,0):C.R2GC_649_193})

V_30 = CTVertex(name = 'V_30',
                type = 'R2',
                particles = [ P.c__tilde__, P.u, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_658_202,(0,1,0):C.R2GC_651_195})

V_31 = CTVertex(name = 'V_31',
                type = 'R2',
                particles = [ P.t__tilde__, P.u, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_783_309,(0,1,0):C.R2GC_745_274})

V_32 = CTVertex(name = 'V_32',
                type = 'R2',
                particles = [ P.t__tilde__, P.u, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_784_310,(0,1,0):C.R2GC_749_278})

V_33 = CTVertex(name = 'V_33',
                type = 'R2',
                particles = [ P.t__tilde__, P.u, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_785_311,(0,1,0):C.R2GC_753_282})

V_34 = CTVertex(name = 'V_34',
                type = 'R2',
                particles = [ P.u__tilde__, P.c, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_653_197,(0,1,0):C.R2GC_648_192})

V_35 = CTVertex(name = 'V_35',
                type = 'R2',
                particles = [ P.u__tilde__, P.c, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_654_198,(0,1,0):C.R2GC_650_194})

V_36 = CTVertex(name = 'V_36',
                type = 'R2',
                particles = [ P.u__tilde__, P.c, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_655_199,(0,1,0):C.R2GC_652_196})

V_37 = CTVertex(name = 'V_37',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_592_146,(0,1,0):C.R2GC_589_143})

V_38 = CTVertex(name = 'V_38',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_593_147,(0,1,0):C.R2GC_590_144})

V_39 = CTVertex(name = 'V_39',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_594_148,(0,1,0):C.R2GC_591_145})

V_40 = CTVertex(name = 'V_40',
                type = 'R2',
                particles = [ P.t__tilde__, P.c, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_786_312,(0,1,0):C.R2GC_746_275})

V_41 = CTVertex(name = 'V_41',
                type = 'R2',
                particles = [ P.t__tilde__, P.c, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_787_313,(0,1,0):C.R2GC_750_279})

V_42 = CTVertex(name = 'V_42',
                type = 'R2',
                particles = [ P.t__tilde__, P.c, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_788_314,(0,1,0):C.R2GC_754_283})

V_43 = CTVertex(name = 'V_43',
                type = 'R2',
                particles = [ P.u__tilde__, P.t, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_775_301,(0,1,0):C.R2GC_747_276})

V_44 = CTVertex(name = 'V_44',
                type = 'R2',
                particles = [ P.u__tilde__, P.t, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_776_302,(0,1,0):C.R2GC_751_280})

V_45 = CTVertex(name = 'V_45',
                type = 'R2',
                particles = [ P.u__tilde__, P.t, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_777_303,(0,1,0):C.R2GC_755_284})

V_46 = CTVertex(name = 'V_46',
                type = 'R2',
                particles = [ P.c__tilde__, P.t, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_780_306,(0,1,0):C.R2GC_748_277})

V_47 = CTVertex(name = 'V_47',
                type = 'R2',
                particles = [ P.c__tilde__, P.t, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_781_307,(0,1,0):C.R2GC_752_281})

V_48 = CTVertex(name = 'V_48',
                type = 'R2',
                particles = [ P.c__tilde__, P.t, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.c, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_782_308,(0,1,0):C.R2GC_756_285})

V_49 = CTVertex(name = 'V_49',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_790_316,(0,1,0):C.R2GC_761_289})

V_50 = CTVertex(name = 'V_50',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_791_317,(0,1,0):C.R2GC_763_290})

V_51 = CTVertex(name = 'V_51',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_792_318,(0,1,0):C.R2GC_765_291})

V_52 = CTVertex(name = 'V_52',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_601_152,(0,1,0):C.R2GC_598_149})

V_53 = CTVertex(name = 'V_53',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_602_153,(0,1,0):C.R2GC_599_150})

V_54 = CTVertex(name = 'V_54',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_603_154,(0,1,0):C.R2GC_600_151})

V_55 = CTVertex(name = 'V_55',
                type = 'R2',
                particles = [ P.s__tilde__, P.d, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_644_188,(0,1,0):C.R2GC_635_179})

V_56 = CTVertex(name = 'V_56',
                type = 'R2',
                particles = [ P.s__tilde__, P.d, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_645_189,(0,1,0):C.R2GC_637_181})

V_57 = CTVertex(name = 'V_57',
                type = 'R2',
                particles = [ P.s__tilde__, P.d, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_646_190,(0,1,0):C.R2GC_639_183})

V_58 = CTVertex(name = 'V_58',
                type = 'R2',
                particles = [ P.b__tilde__, P.d, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_709_244,(0,1,0):C.R2GC_680_218})

V_59 = CTVertex(name = 'V_59',
                type = 'R2',
                particles = [ P.b__tilde__, P.d, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_710_245,(0,1,0):C.R2GC_684_222})

V_60 = CTVertex(name = 'V_60',
                type = 'R2',
                particles = [ P.b__tilde__, P.d, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_711_246,(0,1,0):C.R2GC_688_226})

V_61 = CTVertex(name = 'V_61',
                type = 'R2',
                particles = [ P.d__tilde__, P.s, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_641_185,(0,1,0):C.R2GC_636_180})

V_62 = CTVertex(name = 'V_62',
                type = 'R2',
                particles = [ P.d__tilde__, P.s, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_642_186,(0,1,0):C.R2GC_638_182})

V_63 = CTVertex(name = 'V_63',
                type = 'R2',
                particles = [ P.d__tilde__, P.s, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.d, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_643_187,(0,1,0):C.R2GC_640_184})

V_64 = CTVertex(name = 'V_64',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_610_158,(0,1,0):C.R2GC_607_155})

V_65 = CTVertex(name = 'V_65',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_611_159,(0,1,0):C.R2GC_608_156})

V_66 = CTVertex(name = 'V_66',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_612_160,(0,1,0):C.R2GC_609_157})

V_67 = CTVertex(name = 'V_67',
                type = 'R2',
                particles = [ P.b__tilde__, P.s, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_712_247,(0,1,0):C.R2GC_681_219})

V_68 = CTVertex(name = 'V_68',
                type = 'R2',
                particles = [ P.b__tilde__, P.s, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_713_248,(0,1,0):C.R2GC_685_223})

V_69 = CTVertex(name = 'V_69',
                type = 'R2',
                particles = [ P.b__tilde__, P.s, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_714_249,(0,1,0):C.R2GC_689_227})

V_70 = CTVertex(name = 'V_70',
                type = 'R2',
                particles = [ P.d__tilde__, P.b, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_701_236,(0,1,0):C.R2GC_682_220})

V_71 = CTVertex(name = 'V_71',
                type = 'R2',
                particles = [ P.d__tilde__, P.b, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_702_237,(0,1,0):C.R2GC_686_224})

V_72 = CTVertex(name = 'V_72',
                type = 'R2',
                particles = [ P.d__tilde__, P.b, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_703_238,(0,1,0):C.R2GC_690_228})

V_73 = CTVertex(name = 'V_73',
                type = 'R2',
                particles = [ P.s__tilde__, P.b, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_704_239,(0,1,0):C.R2GC_683_221})

V_74 = CTVertex(name = 'V_74',
                type = 'R2',
                particles = [ P.s__tilde__, P.b, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_705_240,(0,1,0):C.R2GC_687_225})

V_75 = CTVertex(name = 'V_75',
                type = 'R2',
                particles = [ P.s__tilde__, P.b, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_706_241,(0,1,0):C.R2GC_691_229})

V_76 = CTVertex(name = 'V_76',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_715_250,(0,1,0):C.R2GC_694_231})

V_77 = CTVertex(name = 'V_77',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_716_251,(0,1,0):C.R2GC_696_232})

V_78 = CTVertex(name = 'V_78',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS4 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_717_252,(0,1,0):C.R2GC_698_233})

V_79 = CTVertex(name = 'V_79',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_554_36,(0,1,0):C.R2GC_530_1})

V_80 = CTVertex(name = 'V_80',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_554_36,(0,1,0):C.R2GC_530_1})

V_81 = CTVertex(name = 'V_81',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_554_36,(0,1,0):C.R2GC_530_1})

V_82 = CTVertex(name = 'V_82',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_556_38,(0,1,0):C.R2GC_532_3})

V_83 = CTVertex(name = 'V_83',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_556_38,(0,1,0):C.R2GC_532_3})

V_84 = CTVertex(name = 'V_84',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_556_38,(0,1,0):C.R2GC_532_3})

V_85 = CTVertex(name = 'V_85',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_585_141})

V_86 = CTVertex(name = 'V_86',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_585_141})

V_87 = CTVertex(name = 'V_87',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_585_141})

V_88 = CTVertex(name = 'V_88',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_585_141})

V_89 = CTVertex(name = 'V_89',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_585_141})

V_90 = CTVertex(name = 'V_90',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_585_141})

V_91 = CTVertex(name = 'V_91',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_660_204})

V_92 = CTVertex(name = 'V_92',
                type = 'R2',
                particles = [ P.s__tilde__, P.u, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.g, P.s, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_666_210})

V_93 = CTVertex(name = 'V_93',
                type = 'R2',
                particles = [ P.d__tilde__, P.c, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_624_168})

V_94 = CTVertex(name = 'V_94',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_630_174})

V_95 = CTVertex(name = 'V_95',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_740_272})

V_96 = CTVertex(name = 'V_96',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_661_205})

V_97 = CTVertex(name = 'V_97',
                type = 'R2',
                particles = [ P.c__tilde__, P.d, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_626_170})

V_98 = CTVertex(name = 'V_98',
                type = 'R2',
                particles = [ P.u__tilde__, P.s, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.g, P.s, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_667_211})

V_99 = CTVertex(name = 'V_99',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_632_176})

V_100 = CTVertex(name = 'V_100',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_769_295})

V_101 = CTVertex(name = 'V_101',
                 type = 'R2',
                 particles = [ P.u__tilde__, P.u, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_553_35,(0,1,0):C.R2GC_531_2})

V_102 = CTVertex(name = 'V_102',
                 type = 'R2',
                 particles = [ P.c__tilde__, P.c, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_553_35,(0,1,0):C.R2GC_531_2})

V_103 = CTVertex(name = 'V_103',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_553_35,(0,1,0):C.R2GC_531_2})

V_104 = CTVertex(name = 'V_104',
                 type = 'R2',
                 particles = [ P.d__tilde__, P.d, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_555_37,(0,1,0):C.R2GC_533_4})

V_105 = CTVertex(name = 'V_105',
                 type = 'R2',
                 particles = [ P.s__tilde__, P.s, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.R2GC_555_37,(0,1,0):C.R2GC_533_4})

V_106 = CTVertex(name = 'V_106',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_555_37,(0,1,0):C.R2GC_533_4})

V_107 = CTVertex(name = 'V_107',
                 type = 'R2',
                 particles = [ P.u__tilde__, P.u ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_586_142})

V_108 = CTVertex(name = 'V_108',
                 type = 'R2',
                 particles = [ P.c__tilde__, P.c ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_586_142})

V_109 = CTVertex(name = 'V_109',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_741_273,(0,2,0):C.R2GC_741_273,(0,1,0):C.R2GC_586_142,(0,3,0):C.R2GC_586_142})

V_110 = CTVertex(name = 'V_110',
                 type = 'R2',
                 particles = [ P.d__tilde__, P.d ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_586_142})

V_111 = CTVertex(name = 'V_111',
                 type = 'R2',
                 particles = [ P.s__tilde__, P.s ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.R2GC_586_142})

V_112 = CTVertex(name = 'V_112',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_676_217,(0,2,0):C.R2GC_676_217,(0,1,0):C.R2GC_586_142,(0,3,0):C.R2GC_586_142})

V_113 = CTVertex(name = 'V_113',
                 type = 'R2',
                 particles = [ P.g, P.g ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VV1, L.VV2, L.VV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ], [ [P.t] ] ],
                 couplings = {(0,0,2):C.R2GC_721_257,(0,1,0):C.R2GC_538_5,(0,1,3):C.R2GC_538_6,(0,2,1):C.R2GC_720_255,(0,2,2):C.R2GC_720_256})

V_114 = CTVertex(name = 'V_114',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_540_9,(0,0,1):C.R2GC_540_10})

V_115 = CTVertex(name = 'V_115',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_541_11,(0,0,1):C.R2GC_541_12})

V_116 = CTVertex(name = 'V_116',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_542_13,(0,0,1):C.R2GC_542_14})

V_117 = CTVertex(name = 'V_117',
                 type = 'R2',
                 particles = [ P.g, P.g, P.Z, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_550_29,(0,0,1):C.R2GC_550_30,(0,1,0):C.R2GC_550_29,(0,1,1):C.R2GC_550_30,(0,2,0):C.R2GC_550_29,(0,2,1):C.R2GC_550_30})

V_118 = CTVertex(name = 'V_118',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.Z ],
                 color = [ 'Identity(2,3)' ],
                 lorentz = [ L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_551_31,(0,0,1):C.R2GC_551_32,(0,1,0):C.R2GC_551_31,(0,1,1):C.R2GC_551_32,(0,2,0):C.R2GC_551_31,(0,2,1):C.R2GC_551_32})

V_119 = CTVertex(name = 'V_119',
                 type = 'R2',
                 particles = [ P.a, P.a, P.g, P.g ],
                 color = [ 'Identity(3,4)' ],
                 lorentz = [ L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_552_33,(0,0,1):C.R2GC_552_34,(0,1,0):C.R2GC_552_33,(0,1,1):C.R2GC_552_34,(0,2,0):C.R2GC_552_33,(0,2,1):C.R2GC_552_34})

V_120 = CTVertex(name = 'V_120',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__minus__, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b, P.t] ], [ [P.c, P.d] ], [ [P.c, P.s] ], [ [P.d, P.u] ], [ [P.s, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_564_42,(0,0,1):C.R2GC_564_43,(0,0,2):C.R2GC_564_44,(0,0,3):C.R2GC_564_45,(0,0,4):C.R2GC_564_46,(0,1,0):C.R2GC_564_42,(0,1,1):C.R2GC_564_43,(0,1,2):C.R2GC_564_44,(0,1,3):C.R2GC_564_45,(0,1,4):C.R2GC_564_46,(0,2,0):C.R2GC_564_42,(0,2,1):C.R2GC_564_43,(0,2,2):C.R2GC_564_44,(0,2,3):C.R2GC_564_45,(0,2,4):C.R2GC_564_46})

V_121 = CTVertex(name = 'V_121',
                 type = 'R2',
                 particles = [ P.g, P.g, P.g, P.Z ],
                 color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                 lorentz = [ L.VVVV1, L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(1,0,0):C.R2GC_547_23,(1,0,1):C.R2GC_547_24,(0,1,0):C.R2GC_546_21,(0,1,1):C.R2GC_546_22,(0,2,0):C.R2GC_546_21,(0,2,1):C.R2GC_546_22,(0,3,0):C.R2GC_546_21,(0,3,1):C.R2GC_546_22})

V_122 = CTVertex(name = 'V_122',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.g ],
                 color = [ 'd(2,3,4)', 'f(2,3,4)' ],
                 lorentz = [ L.VVVV1, L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(1,0,0):C.R2GC_549_27,(1,0,1):C.R2GC_549_28,(0,1,0):C.R2GC_548_25,(0,1,1):C.R2GC_548_26,(0,2,0):C.R2GC_548_25,(0,2,1):C.R2GC_548_26,(0,3,0):C.R2GC_548_25,(0,3,1):C.R2GC_548_26})

V_123 = CTVertex(name = 'V_123',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G0, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_539_7,(0,0,1):C.R2GC_539_8})

V_124 = CTVertex(name = 'V_124',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G__minus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_559_39})

V_125 = CTVertex(name = 'V_125',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G__minus__, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_560_40})

V_126 = CTVertex(name = 'V_126',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G__plus__, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_561_41})

V_127 = CTVertex(name = 'V_127',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h__minus__, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b, P.c] ], [ [P.b, P.t] ], [ [P.b, P.u] ], [ [P.c, P.d] ], [ [P.c, P.s] ], [ [P.d, P.t] ], [ [P.d, P.u] ], [ [P.s, P.t] ], [ [P.s, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_565_47,(0,0,1):C.R2GC_565_48,(0,0,2):C.R2GC_565_49,(0,0,3):C.R2GC_565_50,(0,0,4):C.R2GC_565_51,(0,0,5):C.R2GC_565_52,(0,0,6):C.R2GC_565_53,(0,0,7):C.R2GC_565_54,(0,0,8):C.R2GC_565_55})

V_128 = CTVertex(name = 'V_128',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G0, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_543_15,(0,0,1):C.R2GC_543_16})

V_129 = CTVertex(name = 'V_129',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h1, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.d] ], [ [P.b, P.s] ], [ [P.c] ], [ [P.c, P.t] ], [ [P.c, P.u] ], [ [P.d] ], [ [P.d, P.s] ], [ [P.s] ], [ [P.t] ], [ [P.t, P.u] ], [ [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_566_56,(0,0,3):C.R2GC_566_57,(0,0,6):C.R2GC_566_58,(0,0,8):C.R2GC_566_59,(0,0,9):C.R2GC_566_60,(0,0,11):C.R2GC_566_61,(0,0,1):C.R2GC_566_62,(0,0,2):C.R2GC_566_63,(0,0,4):C.R2GC_566_64,(0,0,5):C.R2GC_566_65,(0,0,7):C.R2GC_566_66,(0,0,10):C.R2GC_566_67})

V_130 = CTVertex(name = 'V_130',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G0, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_544_17,(0,0,1):C.R2GC_544_18})

V_131 = CTVertex(name = 'V_131',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h1, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.d] ], [ [P.b, P.s] ], [ [P.c] ], [ [P.c, P.t] ], [ [P.c, P.u] ], [ [P.d] ], [ [P.d, P.s] ], [ [P.s] ], [ [P.t] ], [ [P.t, P.u] ], [ [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_567_68,(0,0,3):C.R2GC_567_69,(0,0,6):C.R2GC_567_70,(0,0,8):C.R2GC_567_71,(0,0,9):C.R2GC_567_72,(0,0,11):C.R2GC_567_73,(0,0,1):C.R2GC_567_74,(0,0,2):C.R2GC_567_75,(0,0,4):C.R2GC_567_76,(0,0,5):C.R2GC_567_77,(0,0,7):C.R2GC_567_78,(0,0,10):C.R2GC_567_79})

V_132 = CTVertex(name = 'V_132',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h2, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.d] ], [ [P.b, P.s] ], [ [P.c] ], [ [P.c, P.t] ], [ [P.c, P.u] ], [ [P.d] ], [ [P.d, P.s] ], [ [P.s] ], [ [P.t] ], [ [P.t, P.u] ], [ [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_568_80,(0,0,3):C.R2GC_568_81,(0,0,6):C.R2GC_568_82,(0,0,8):C.R2GC_568_83,(0,0,9):C.R2GC_568_84,(0,0,11):C.R2GC_568_85,(0,0,1):C.R2GC_568_86,(0,0,2):C.R2GC_568_87,(0,0,4):C.R2GC_568_88,(0,0,5):C.R2GC_568_89,(0,0,7):C.R2GC_568_90,(0,0,10):C.R2GC_568_91})

V_133 = CTVertex(name = 'V_133',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G0, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_545_19,(0,0,1):C.R2GC_545_20})

V_134 = CTVertex(name = 'V_134',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h1, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.d] ], [ [P.b, P.s] ], [ [P.c] ], [ [P.c, P.t] ], [ [P.c, P.u] ], [ [P.d] ], [ [P.d, P.s] ], [ [P.s] ], [ [P.t] ], [ [P.t, P.u] ], [ [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_569_92,(0,0,3):C.R2GC_569_93,(0,0,6):C.R2GC_569_94,(0,0,8):C.R2GC_569_95,(0,0,9):C.R2GC_569_96,(0,0,11):C.R2GC_569_97,(0,0,1):C.R2GC_569_98,(0,0,2):C.R2GC_569_99,(0,0,4):C.R2GC_569_100,(0,0,5):C.R2GC_569_101,(0,0,7):C.R2GC_569_102,(0,0,10):C.R2GC_569_103})

V_135 = CTVertex(name = 'V_135',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h2, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.d] ], [ [P.b, P.s] ], [ [P.c] ], [ [P.c, P.t] ], [ [P.c, P.u] ], [ [P.d] ], [ [P.d, P.s] ], [ [P.s] ], [ [P.t] ], [ [P.t, P.u] ], [ [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_570_104,(0,0,3):C.R2GC_570_105,(0,0,6):C.R2GC_570_106,(0,0,8):C.R2GC_570_107,(0,0,9):C.R2GC_570_108,(0,0,11):C.R2GC_570_109,(0,0,1):C.R2GC_570_110,(0,0,2):C.R2GC_570_111,(0,0,4):C.R2GC_570_112,(0,0,5):C.R2GC_570_113,(0,0,7):C.R2GC_570_114,(0,0,10):C.R2GC_570_115})

V_136 = CTVertex(name = 'V_136',
                 type = 'R2',
                 particles = [ P.g, P.g, P.h3, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.d] ], [ [P.b, P.s] ], [ [P.c] ], [ [P.c, P.t] ], [ [P.c, P.u] ], [ [P.d] ], [ [P.d, P.s] ], [ [P.s] ], [ [P.t] ], [ [P.t, P.u] ], [ [P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_571_116,(0,0,3):C.R2GC_571_117,(0,0,6):C.R2GC_571_118,(0,0,8):C.R2GC_571_119,(0,0,9):C.R2GC_571_120,(0,0,11):C.R2GC_571_121,(0,0,1):C.R2GC_571_122,(0,0,2):C.R2GC_571_123,(0,0,4):C.R2GC_571_124,(0,0,5):C.R2GC_571_125,(0,0,7):C.R2GC_571_126,(0,0,10):C.R2GC_571_127})

V_137 = CTVertex(name = 'V_137',
                 type = 'UV',
                 particles = [ P.g, P.g, P.g ],
                 color = [ 'f(1,2,3)' ],
                 lorentz = [ L.VVV1, L.VVV2, L.VVV3, L.VVV4, L.VVV5, L.VVV6, L.VVV7, L.VVV8 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_722_265,(0,0,1):C.UVGC_722_266,(0,0,2):C.UVGC_722_267,(0,0,3):C.UVGC_577_5,(0,0,4):C.UVGC_722_268,(0,1,0):C.UVGC_725_273,(0,1,1):C.UVGC_725_274,(0,1,2):C.UVGC_725_275,(0,1,3):C.UVGC_725_276,(0,1,4):C.UVGC_725_277,(0,3,0):C.UVGC_725_273,(0,3,1):C.UVGC_725_274,(0,3,2):C.UVGC_727_280,(0,3,3):C.UVGC_576_3,(0,3,4):C.UVGC_725_277,(0,5,0):C.UVGC_722_265,(0,5,1):C.UVGC_722_266,(0,5,2):C.UVGC_724_271,(0,5,3):C.UVGC_724_272,(0,5,4):C.UVGC_722_268,(0,6,0):C.UVGC_722_265,(0,6,1):C.UVGC_722_266,(0,6,2):C.UVGC_723_269,(0,6,3):C.UVGC_723_270,(0,6,4):C.UVGC_722_268,(0,7,0):C.UVGC_725_273,(0,7,1):C.UVGC_725_274,(0,7,2):C.UVGC_726_278,(0,7,3):C.UVGC_726_279,(0,7,4):C.UVGC_725_277,(0,2,2):C.UVGC_576_2,(0,2,3):C.UVGC_576_3,(0,4,2):C.UVGC_577_4,(0,4,3):C.UVGC_577_5})

V_138 = CTVertex(name = 'V_138',
                 type = 'UV',
                 particles = [ P.g, P.g, P.g, P.g ],
                 color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
                 lorentz = [ L.VVVV2, L.VVVV3, L.VVVV4 ],
                 loop_particles = [ [ [P.b] ], [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(2,0,3):C.UVGC_579_9,(2,0,4):C.UVGC_579_8,(0,0,3):C.UVGC_579_9,(0,0,4):C.UVGC_579_8,(4,0,3):C.UVGC_578_6,(4,0,4):C.UVGC_578_7,(3,0,3):C.UVGC_578_6,(3,0,4):C.UVGC_578_7,(8,0,3):C.UVGC_579_8,(8,0,4):C.UVGC_579_9,(7,0,0):C.UVGC_731_292,(7,0,2):C.UVGC_731_293,(7,0,3):C.UVGC_731_294,(7,0,4):C.UVGC_731_295,(7,0,5):C.UVGC_731_296,(6,0,0):C.UVGC_731_292,(6,0,2):C.UVGC_731_293,(6,0,3):C.UVGC_732_297,(6,0,4):C.UVGC_732_298,(6,0,5):C.UVGC_731_296,(5,0,3):C.UVGC_578_6,(5,0,4):C.UVGC_578_7,(1,0,3):C.UVGC_578_6,(1,0,4):C.UVGC_578_7,(11,0,3):C.UVGC_582_12,(11,0,4):C.UVGC_582_13,(10,0,3):C.UVGC_582_12,(10,0,4):C.UVGC_582_13,(9,0,3):C.UVGC_581_10,(9,0,4):C.UVGC_581_11,(2,1,3):C.UVGC_579_9,(2,1,4):C.UVGC_579_8,(0,1,3):C.UVGC_579_9,(0,1,4):C.UVGC_579_8,(4,1,3):C.UVGC_578_6,(4,1,4):C.UVGC_578_7,(3,1,3):C.UVGC_578_6,(3,1,4):C.UVGC_578_7,(8,1,0):C.UVGC_733_299,(8,1,2):C.UVGC_733_300,(8,1,3):C.UVGC_733_301,(8,1,4):C.UVGC_733_302,(8,1,5):C.UVGC_733_303,(6,1,0):C.UVGC_728_281,(6,1,3):C.UVGC_728_282,(6,1,4):C.UVGC_728_283,(6,1,5):C.UVGC_728_284,(7,1,1):C.UVGC_583_14,(7,1,3):C.UVGC_584_16,(7,1,4):C.UVGC_584_17,(5,1,3):C.UVGC_578_6,(5,1,4):C.UVGC_578_7,(1,1,3):C.UVGC_578_6,(1,1,4):C.UVGC_578_7,(11,1,3):C.UVGC_582_12,(11,1,4):C.UVGC_582_13,(10,1,3):C.UVGC_582_12,(10,1,4):C.UVGC_582_13,(9,1,3):C.UVGC_581_10,(9,1,4):C.UVGC_581_11,(2,2,3):C.UVGC_579_9,(2,2,4):C.UVGC_579_8,(0,2,3):C.UVGC_579_9,(0,2,4):C.UVGC_579_8,(4,2,3):C.UVGC_578_6,(4,2,4):C.UVGC_578_7,(3,2,3):C.UVGC_578_6,(3,2,4):C.UVGC_578_7,(8,2,0):C.UVGC_730_287,(8,2,2):C.UVGC_730_288,(8,2,3):C.UVGC_730_289,(8,2,4):C.UVGC_730_290,(8,2,5):C.UVGC_730_291,(6,2,1):C.UVGC_583_14,(6,2,3):C.UVGC_583_15,(6,2,4):C.UVGC_581_10,(7,2,0):C.UVGC_728_281,(7,2,3):C.UVGC_729_285,(7,2,4):C.UVGC_729_286,(7,2,5):C.UVGC_728_284,(5,2,3):C.UVGC_578_6,(5,2,4):C.UVGC_578_7,(1,2,3):C.UVGC_578_6,(1,2,4):C.UVGC_578_7,(11,2,3):C.UVGC_582_12,(11,2,4):C.UVGC_582_13,(10,2,3):C.UVGC_582_12,(10,2,4):C.UVGC_582_13,(9,2,3):C.UVGC_581_10,(9,2,4):C.UVGC_581_11})

V_139 = CTVertex(name = 'V_139',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_664_132,(0,0,1):C.UVGC_664_133,(0,1,0):C.UVGC_659_122,(0,1,1):C.UVGC_659_123})

V_140 = CTVertex(name = 'V_140',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.u, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s], [P.g, P.u] ], [ [P.g, P.s, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_670_144,(0,0,1):C.UVGC_670_145,(0,1,0):C.UVGC_665_134,(0,1,1):C.UVGC_665_135})

V_141 = CTVertex(name = 'V_141',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.u, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_718_253,(0,0,2):C.UVGC_718_254,(0,0,1):C.UVGC_718_255,(0,1,0):C.UVGC_674_149,(0,1,2):C.UVGC_674_150,(0,1,1):C.UVGC_674_151})

V_142 = CTVertex(name = 'V_142',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.c, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.d, P.g] ], [ [P.c, P.g], [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_628_60,(0,0,0):C.UVGC_628_61,(0,1,1):C.UVGC_623_50,(0,1,0):C.UVGC_623_51})

V_143 = CTVertex(name = 'V_143',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_634_72,(0,0,1):C.UVGC_634_73,(0,1,0):C.UVGC_629_62,(0,1,1):C.UVGC_629_63})

V_144 = CTVertex(name = 'V_144',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.c, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.c, P.g] ], [ [P.b, P.g] ], [ [P.c, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_719_256,(0,0,2):C.UVGC_719_257,(0,0,0):C.UVGC_719_258,(0,1,1):C.UVGC_675_152,(0,1,2):C.UVGC_675_153,(0,1,0):C.UVGC_675_154})

V_145 = CTVertex(name = 'V_145',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.t, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g] ], [ [P.d, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_778_408,(0,0,2):C.UVGC_778_409,(0,0,1):C.UVGC_778_410,(0,1,0):C.UVGC_737_307,(0,1,2):C.UVGC_737_308,(0,1,1):C.UVGC_737_309})

V_146 = CTVertex(name = 'V_146',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.t, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s] ], [ [P.g, P.s, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_779_411,(0,0,2):C.UVGC_779_412,(0,0,1):C.UVGC_779_413,(0,1,0):C.UVGC_738_310,(0,1,2):C.UVGC_738_311,(0,1,1):C.UVGC_738_312})

V_147 = CTVertex(name = 'V_147',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.h__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_789_441,(0,0,2):C.UVGC_789_442,(0,0,1):C.UVGC_789_443,(0,1,0):C.UVGC_739_313,(0,1,2):C.UVGC_739_314,(0,1,1):C.UVGC_739_315})

V_148 = CTVertex(name = 'V_148',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_663_130,(0,0,1):C.UVGC_663_131,(0,1,0):C.UVGC_662_128,(0,1,1):C.UVGC_662_129})

V_149 = CTVertex(name = 'V_149',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.d, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.d, P.g] ], [ [P.c, P.g], [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_627_58,(0,0,0):C.UVGC_627_59,(0,1,1):C.UVGC_625_54,(0,1,0):C.UVGC_625_55})

V_150 = CTVertex(name = 'V_150',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.d, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g] ], [ [P.d, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_772_390,(0,0,2):C.UVGC_772_391,(0,0,1):C.UVGC_772_392,(0,1,0):C.UVGC_766_372,(0,1,2):C.UVGC_766_373,(0,1,1):C.UVGC_766_374})

V_151 = CTVertex(name = 'V_151',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.s, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s], [P.g, P.u] ], [ [P.g, P.s, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_669_142,(0,0,1):C.UVGC_669_143,(0,1,0):C.UVGC_668_140,(0,1,1):C.UVGC_668_141})

V_152 = CTVertex(name = 'V_152',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_633_70,(0,0,1):C.UVGC_633_71,(0,1,0):C.UVGC_631_66,(0,1,1):C.UVGC_631_67})

V_153 = CTVertex(name = 'V_153',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.s, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s] ], [ [P.g, P.s, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_773_393,(0,0,2):C.UVGC_773_394,(0,0,1):C.UVGC_773_395,(0,1,0):C.UVGC_767_375,(0,1,2):C.UVGC_767_376,(0,1,1):C.UVGC_767_377})

V_154 = CTVertex(name = 'V_154',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.b, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_707_226,(0,0,2):C.UVGC_707_227,(0,0,1):C.UVGC_707_228,(0,1,0):C.UVGC_699_202,(0,1,2):C.UVGC_699_203,(0,1,1):C.UVGC_699_204})

V_155 = CTVertex(name = 'V_155',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.b, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.c, P.g] ], [ [P.b, P.g] ], [ [P.c, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_708_229,(0,0,2):C.UVGC_708_230,(0,0,0):C.UVGC_708_231,(0,1,1):C.UVGC_700_205,(0,1,2):C.UVGC_700_206,(0,1,0):C.UVGC_700_207})

V_156 = CTVertex(name = 'V_156',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.h__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_774_396,(0,0,2):C.UVGC_774_397,(0,0,1):C.UVGC_774_398,(0,1,0):C.UVGC_768_378,(0,1,2):C.UVGC_768_379,(0,1,1):C.UVGC_768_380})

V_157 = CTVertex(name = 'V_157',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_770_384,(0,0,2):C.UVGC_770_385,(0,0,1):C.UVGC_770_386,(0,1,0):C.UVGC_771_387,(0,1,2):C.UVGC_771_388,(0,1,1):C.UVGC_771_389})

V_158 = CTVertex(name = 'V_158',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS1 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_692_195})

V_159 = CTVertex(name = 'V_159',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_759_363,(0,0,2):C.UVGC_759_364,(0,0,1):C.UVGC_759_365,(0,1,0):C.UVGC_757_359,(0,1,2):C.UVGC_757_360,(0,1,1):C.UVGC_757_361})

V_160 = CTVertex(name = 'V_160',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS1 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_758_362})

V_161 = CTVertex(name = 'V_161',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_620_47,(0,1,0):C.UVGC_617_44})

V_162 = CTVertex(name = 'V_162',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_621_48,(0,1,0):C.UVGC_618_45})

V_163 = CTVertex(name = 'V_163',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_622_49,(0,1,0):C.UVGC_619_46})

V_164 = CTVertex(name = 'V_164',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.u, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.u] ], [ [P.c, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_656_116,(0,0,1):C.UVGC_656_117,(0,1,0):C.UVGC_647_98,(0,1,1):C.UVGC_647_99})

V_165 = CTVertex(name = 'V_165',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.u, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.u] ], [ [P.c, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_657_118,(0,0,1):C.UVGC_657_119,(0,1,0):C.UVGC_649_102,(0,1,1):C.UVGC_649_103})

V_166 = CTVertex(name = 'V_166',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.u, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.u] ], [ [P.c, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_658_120,(0,0,1):C.UVGC_658_121,(0,1,0):C.UVGC_651_106,(0,1,1):C.UVGC_651_107})

V_167 = CTVertex(name = 'V_167',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.u, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ], [ [P.g, P.t, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_783_423,(0,0,2):C.UVGC_783_424,(0,0,1):C.UVGC_783_425,(0,1,0):C.UVGC_745_323,(0,1,2):C.UVGC_745_324,(0,1,1):C.UVGC_745_325})

V_168 = CTVertex(name = 'V_168',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.u, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ], [ [P.g, P.t, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_784_426,(0,0,2):C.UVGC_784_427,(0,0,1):C.UVGC_784_428,(0,1,0):C.UVGC_749_335,(0,1,2):C.UVGC_749_336,(0,1,1):C.UVGC_749_337})

V_169 = CTVertex(name = 'V_169',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.u, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ], [ [P.g, P.t, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_785_429,(0,0,2):C.UVGC_785_430,(0,0,1):C.UVGC_785_431,(0,1,0):C.UVGC_753_347,(0,1,2):C.UVGC_753_348,(0,1,1):C.UVGC_753_349})

V_170 = CTVertex(name = 'V_170',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.c, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.u] ], [ [P.c, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_653_110,(0,0,1):C.UVGC_653_111,(0,1,0):C.UVGC_648_100,(0,1,1):C.UVGC_648_101})

V_171 = CTVertex(name = 'V_171',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.c, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.u] ], [ [P.c, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_654_112,(0,0,1):C.UVGC_654_113,(0,1,0):C.UVGC_650_104,(0,1,1):C.UVGC_650_105})

V_172 = CTVertex(name = 'V_172',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.c, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.u] ], [ [P.c, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_655_114,(0,0,1):C.UVGC_655_115,(0,1,0):C.UVGC_652_108,(0,1,1):C.UVGC_652_109})

V_173 = CTVertex(name = 'V_173',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_592_29,(0,1,0):C.UVGC_589_26})

V_174 = CTVertex(name = 'V_174',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_593_30,(0,1,0):C.UVGC_590_27})

V_175 = CTVertex(name = 'V_175',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_594_31,(0,1,0):C.UVGC_591_28})

V_176 = CTVertex(name = 'V_176',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.c, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ], [ [P.c, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_786_432,(0,0,2):C.UVGC_786_433,(0,0,1):C.UVGC_786_434,(0,1,0):C.UVGC_746_326,(0,1,2):C.UVGC_746_327,(0,1,1):C.UVGC_746_328})

V_177 = CTVertex(name = 'V_177',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.c, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ], [ [P.c, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_787_435,(0,0,2):C.UVGC_787_436,(0,0,1):C.UVGC_787_437,(0,1,0):C.UVGC_750_338,(0,1,2):C.UVGC_750_339,(0,1,1):C.UVGC_750_340})

V_178 = CTVertex(name = 'V_178',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.c, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ], [ [P.c, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_788_438,(0,0,2):C.UVGC_788_439,(0,0,1):C.UVGC_788_440,(0,1,0):C.UVGC_754_350,(0,1,2):C.UVGC_754_351,(0,1,1):C.UVGC_754_352})

V_179 = CTVertex(name = 'V_179',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.t, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ], [ [P.g, P.t, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_775_399,(0,0,2):C.UVGC_775_400,(0,0,1):C.UVGC_775_401,(0,1,0):C.UVGC_747_329,(0,1,2):C.UVGC_747_330,(0,1,1):C.UVGC_747_331})

V_180 = CTVertex(name = 'V_180',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.t, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ], [ [P.g, P.t, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_776_402,(0,0,2):C.UVGC_776_403,(0,0,1):C.UVGC_776_404,(0,1,0):C.UVGC_751_341,(0,1,2):C.UVGC_751_342,(0,1,1):C.UVGC_751_343})

V_181 = CTVertex(name = 'V_181',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.t, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ], [ [P.g, P.t, P.u] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_777_405,(0,0,2):C.UVGC_777_406,(0,0,1):C.UVGC_777_407,(0,1,0):C.UVGC_755_353,(0,1,2):C.UVGC_755_354,(0,1,1):C.UVGC_755_355})

V_182 = CTVertex(name = 'V_182',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.t, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ], [ [P.c, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_780_414,(0,0,2):C.UVGC_780_415,(0,0,1):C.UVGC_780_416,(0,1,0):C.UVGC_748_332,(0,1,2):C.UVGC_748_333,(0,1,1):C.UVGC_748_334})

V_183 = CTVertex(name = 'V_183',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.t, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ], [ [P.c, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_781_417,(0,0,2):C.UVGC_781_418,(0,0,1):C.UVGC_781_419,(0,1,0):C.UVGC_752_344,(0,1,2):C.UVGC_752_345,(0,1,1):C.UVGC_752_346})

V_184 = CTVertex(name = 'V_184',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.t, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.c, P.g] ], [ [P.c, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_782_420,(0,0,2):C.UVGC_782_421,(0,0,1):C.UVGC_782_422,(0,1,0):C.UVGC_756_356,(0,1,2):C.UVGC_756_357,(0,1,1):C.UVGC_756_358})

V_185 = CTVertex(name = 'V_185',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2, L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_790_444,(0,2,0):C.UVGC_761_367,(0,0,0):C.UVGC_760_366})

V_186 = CTVertex(name = 'V_186',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2, L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_791_445,(0,2,0):C.UVGC_763_369,(0,0,0):C.UVGC_762_368})

V_187 = CTVertex(name = 'V_187',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2, L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_792_446,(0,2,0):C.UVGC_765_371,(0,0,0):C.UVGC_764_370})

V_188 = CTVertex(name = 'V_188',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_601_35,(0,1,0):C.UVGC_598_32})

V_189 = CTVertex(name = 'V_189',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_602_36,(0,1,0):C.UVGC_599_33})

V_190 = CTVertex(name = 'V_190',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_603_37,(0,1,0):C.UVGC_600_34})

V_191 = CTVertex(name = 'V_191',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.d, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.s] ], [ [P.d, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_644_92,(0,0,1):C.UVGC_644_93,(0,1,0):C.UVGC_635_74,(0,1,1):C.UVGC_635_75})

V_192 = CTVertex(name = 'V_192',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.d, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.s] ], [ [P.d, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_645_94,(0,0,1):C.UVGC_645_95,(0,1,0):C.UVGC_637_78,(0,1,1):C.UVGC_637_79})

V_193 = CTVertex(name = 'V_193',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.d, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.s] ], [ [P.d, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_646_96,(0,0,1):C.UVGC_646_97,(0,1,0):C.UVGC_639_82,(0,1,1):C.UVGC_639_83})

V_194 = CTVertex(name = 'V_194',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.d, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.d, P.g] ], [ [P.b, P.g] ], [ [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_709_232,(0,0,2):C.UVGC_709_233,(0,0,0):C.UVGC_709_234,(0,1,1):C.UVGC_680_159,(0,1,2):C.UVGC_680_160,(0,1,0):C.UVGC_680_161})

V_195 = CTVertex(name = 'V_195',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.d, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.d, P.g] ], [ [P.b, P.g] ], [ [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_710_235,(0,0,2):C.UVGC_710_236,(0,0,0):C.UVGC_710_237,(0,1,1):C.UVGC_684_171,(0,1,2):C.UVGC_684_172,(0,1,0):C.UVGC_684_173})

V_196 = CTVertex(name = 'V_196',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.d, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.d, P.g] ], [ [P.b, P.g] ], [ [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_711_238,(0,0,2):C.UVGC_711_239,(0,0,0):C.UVGC_711_240,(0,1,1):C.UVGC_688_183,(0,1,2):C.UVGC_688_184,(0,1,0):C.UVGC_688_185})

V_197 = CTVertex(name = 'V_197',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.s, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.s] ], [ [P.d, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_641_86,(0,0,1):C.UVGC_641_87,(0,1,0):C.UVGC_636_76,(0,1,1):C.UVGC_636_77})

V_198 = CTVertex(name = 'V_198',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.s, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.s] ], [ [P.d, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_642_88,(0,0,1):C.UVGC_642_89,(0,1,0):C.UVGC_638_80,(0,1,1):C.UVGC_638_81})

V_199 = CTVertex(name = 'V_199',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.s, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.s] ], [ [P.d, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_643_90,(0,0,1):C.UVGC_643_91,(0,1,0):C.UVGC_640_84,(0,1,1):C.UVGC_640_85})

V_200 = CTVertex(name = 'V_200',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_610_41,(0,1,0):C.UVGC_607_38})

V_201 = CTVertex(name = 'V_201',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_611_42,(0,1,0):C.UVGC_608_39})

V_202 = CTVertex(name = 'V_202',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_612_43,(0,1,0):C.UVGC_609_40})

V_203 = CTVertex(name = 'V_203',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.s, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.s] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_712_241,(0,0,2):C.UVGC_712_242,(0,0,1):C.UVGC_712_243,(0,1,0):C.UVGC_681_162,(0,1,2):C.UVGC_681_163,(0,1,1):C.UVGC_681_164})

V_204 = CTVertex(name = 'V_204',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.s, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.s] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_713_244,(0,0,2):C.UVGC_713_245,(0,0,1):C.UVGC_713_246,(0,1,0):C.UVGC_685_174,(0,1,2):C.UVGC_685_175,(0,1,1):C.UVGC_685_176})

V_205 = CTVertex(name = 'V_205',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.s, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.s] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_714_247,(0,0,2):C.UVGC_714_248,(0,0,1):C.UVGC_714_249,(0,1,0):C.UVGC_689_186,(0,1,2):C.UVGC_689_187,(0,1,1):C.UVGC_689_188})

V_206 = CTVertex(name = 'V_206',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.b, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.d, P.g] ], [ [P.b, P.g] ], [ [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_701_208,(0,0,2):C.UVGC_701_209,(0,0,0):C.UVGC_701_210,(0,1,1):C.UVGC_682_165,(0,1,2):C.UVGC_682_166,(0,1,0):C.UVGC_682_167})

V_207 = CTVertex(name = 'V_207',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.b, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.d, P.g] ], [ [P.b, P.g] ], [ [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_702_211,(0,0,2):C.UVGC_702_212,(0,0,0):C.UVGC_702_213,(0,1,1):C.UVGC_686_177,(0,1,2):C.UVGC_686_178,(0,1,0):C.UVGC_686_179})

V_208 = CTVertex(name = 'V_208',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.b, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.d, P.g] ], [ [P.b, P.g] ], [ [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_703_214,(0,0,2):C.UVGC_703_215,(0,0,0):C.UVGC_703_216,(0,1,1):C.UVGC_690_189,(0,1,2):C.UVGC_690_190,(0,1,0):C.UVGC_690_191})

V_209 = CTVertex(name = 'V_209',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.b, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.s] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_704_217,(0,0,2):C.UVGC_704_218,(0,0,1):C.UVGC_704_219,(0,1,0):C.UVGC_683_168,(0,1,2):C.UVGC_683_169,(0,1,1):C.UVGC_683_170})

V_210 = CTVertex(name = 'V_210',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.b, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.s] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_705_220,(0,0,2):C.UVGC_705_221,(0,0,1):C.UVGC_705_222,(0,1,0):C.UVGC_687_180,(0,1,2):C.UVGC_687_181,(0,1,1):C.UVGC_687_182})

V_211 = CTVertex(name = 'V_211',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.b, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.s] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_706_223,(0,0,2):C.UVGC_706_224,(0,0,1):C.UVGC_706_225,(0,1,0):C.UVGC_691_192,(0,1,2):C.UVGC_691_193,(0,1,1):C.UVGC_691_194})

V_212 = CTVertex(name = 'V_212',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.h1 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2, L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,1,0):C.UVGC_715_250,(0,2,0):C.UVGC_694_197,(0,0,0):C.UVGC_693_196})

V_213 = CTVertex(name = 'V_213',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.h2 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2, L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,1,0):C.UVGC_716_251,(0,2,0):C.UVGC_696_199,(0,0,0):C.UVGC_695_198})

V_214 = CTVertex(name = 'V_214',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.h3 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2, L.FFS3, L.FFS4 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,1,0):C.UVGC_717_252,(0,2,0):C.UVGC_698_201,(0,0,0):C.UVGC_697_200})

V_215 = CTVertex(name = 'V_215',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_744_322,(0,1,0):C.UVGC_736_306})

V_216 = CTVertex(name = 'V_216',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.a ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_679_158,(0,1,0):C.UVGC_673_148})

V_217 = CTVertex(name = 'V_217',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.u] ], [ [P.t] ] ],
                 couplings = {(0,0,4):C.UVGC_585_18,(0,1,0):C.UVGC_587_20,(0,1,1):C.UVGC_587_21,(0,1,2):C.UVGC_587_22,(0,1,3):C.UVGC_587_23,(0,1,5):C.UVGC_587_24,(0,1,4):C.UVGC_587_25,(0,2,0):C.UVGC_587_20,(0,2,1):C.UVGC_587_21,(0,2,2):C.UVGC_587_22,(0,2,3):C.UVGC_587_23,(0,2,5):C.UVGC_587_24,(0,2,4):C.UVGC_587_25})

V_218 = CTVertex(name = 'V_218',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.c, P.g] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,2):C.UVGC_585_18,(0,1,0):C.UVGC_587_20,(0,1,1):C.UVGC_587_21,(0,1,3):C.UVGC_587_22,(0,1,4):C.UVGC_587_23,(0,1,5):C.UVGC_587_24,(0,1,2):C.UVGC_587_25,(0,2,0):C.UVGC_587_20,(0,2,1):C.UVGC_587_21,(0,2,3):C.UVGC_587_22,(0,2,4):C.UVGC_587_23,(0,2,5):C.UVGC_587_24,(0,2,2):C.UVGC_587_25})

V_219 = CTVertex(name = 'V_219',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ], [ [P.t] ] ],
                 couplings = {(0,0,4):C.UVGC_585_18,(0,1,0):C.UVGC_587_20,(0,1,1):C.UVGC_587_21,(0,1,2):C.UVGC_587_22,(0,1,3):C.UVGC_587_23,(0,1,5):C.UVGC_587_24,(0,1,4):C.UVGC_735_305,(0,2,0):C.UVGC_587_20,(0,2,1):C.UVGC_587_21,(0,2,2):C.UVGC_587_22,(0,2,3):C.UVGC_587_23,(0,2,5):C.UVGC_587_24,(0,2,4):C.UVGC_735_305})

V_220 = CTVertex(name = 'V_220',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.d, P.g] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,2):C.UVGC_585_18,(0,1,0):C.UVGC_587_20,(0,1,1):C.UVGC_587_21,(0,1,3):C.UVGC_587_22,(0,1,4):C.UVGC_587_23,(0,1,5):C.UVGC_587_24,(0,1,2):C.UVGC_587_25,(0,2,0):C.UVGC_587_20,(0,2,1):C.UVGC_587_21,(0,2,3):C.UVGC_587_22,(0,2,4):C.UVGC_587_23,(0,2,5):C.UVGC_587_24,(0,2,2):C.UVGC_587_25})

V_221 = CTVertex(name = 'V_221',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,4):C.UVGC_585_18,(0,1,0):C.UVGC_587_20,(0,1,1):C.UVGC_587_21,(0,1,2):C.UVGC_587_22,(0,1,3):C.UVGC_587_23,(0,1,5):C.UVGC_587_24,(0,1,4):C.UVGC_587_25,(0,2,0):C.UVGC_587_20,(0,2,1):C.UVGC_587_21,(0,2,2):C.UVGC_587_22,(0,2,3):C.UVGC_587_23,(0,2,5):C.UVGC_587_24,(0,2,4):C.UVGC_587_25})

V_222 = CTVertex(name = 'V_222',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.g] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,1):C.UVGC_585_18,(0,1,0):C.UVGC_587_20,(0,1,2):C.UVGC_587_21,(0,1,3):C.UVGC_587_22,(0,1,4):C.UVGC_587_23,(0,1,5):C.UVGC_587_24,(0,1,1):C.UVGC_672_147,(0,2,0):C.UVGC_587_20,(0,2,2):C.UVGC_587_21,(0,2,3):C.UVGC_587_22,(0,2,4):C.UVGC_587_23,(0,2,5):C.UVGC_587_24,(0,2,1):C.UVGC_672_147})

V_223 = CTVertex(name = 'V_223',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_660_124,(0,0,1):C.UVGC_660_125})

V_224 = CTVertex(name = 'V_224',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.u, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.g, P.s], [P.g, P.u] ], [ [P.g, P.s, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_666_136,(0,0,1):C.UVGC_666_137})

V_225 = CTVertex(name = 'V_225',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.c, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.c, P.d, P.g] ], [ [P.c, P.g], [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_624_52,(0,0,0):C.UVGC_624_53})

V_226 = CTVertex(name = 'V_226',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_630_64,(0,0,1):C.UVGC_630_65})

V_227 = CTVertex(name = 'V_227',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_740_316,(0,0,2):C.UVGC_740_317,(0,0,1):C.UVGC_740_318})

V_228 = CTVertex(name = 'V_228',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_661_126,(0,0,1):C.UVGC_661_127})

V_229 = CTVertex(name = 'V_229',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.d, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.c, P.d, P.g] ], [ [P.c, P.g], [P.d, P.g] ] ],
                 couplings = {(0,0,1):C.UVGC_626_56,(0,0,0):C.UVGC_626_57})

V_230 = CTVertex(name = 'V_230',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.s, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.g, P.s], [P.g, P.u] ], [ [P.g, P.s, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_667_138,(0,0,1):C.UVGC_667_139})

V_231 = CTVertex(name = 'V_231',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_632_68,(0,0,1):C.UVGC_632_69})

V_232 = CTVertex(name = 'V_232',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_769_381,(0,0,2):C.UVGC_769_382,(0,0,1):C.UVGC_769_383})

V_233 = CTVertex(name = 'V_233',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_742_320,(0,1,0):C.UVGC_743_321})

V_234 = CTVertex(name = 'V_234',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV3 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_677_156,(0,1,0):C.UVGC_678_157})

V_235 = CTVertex(name = 'V_235',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1, L.FF3, L.FF5 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_586_19,(0,1,0):C.UVGC_572_1,(0,2,0):C.UVGC_572_1})

V_236 = CTVertex(name = 'V_236',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1, L.FF3, L.FF5 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_586_19,(0,1,0):C.UVGC_572_1,(0,2,0):C.UVGC_572_1})

V_237 = CTVertex(name = 'V_237',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_741_319,(0,2,0):C.UVGC_741_319,(0,1,0):C.UVGC_734_304,(0,3,0):C.UVGC_734_304})

V_238 = CTVertex(name = 'V_238',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1, L.FF3, L.FF5 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_586_19,(0,1,0):C.UVGC_572_1,(0,2,0):C.UVGC_572_1})

V_239 = CTVertex(name = 'V_239',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1, L.FF3, L.FF5 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_586_19,(0,1,0):C.UVGC_572_1,(0,2,0):C.UVGC_572_1})

V_240 = CTVertex(name = 'V_240',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_676_155,(0,2,0):C.UVGC_676_155,(0,1,0):C.UVGC_671_146,(0,3,0):C.UVGC_671_146})

V_241 = CTVertex(name = 'V_241',
                 type = 'UV',
                 particles = [ P.g, P.g ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VV1, L.VV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_721_261,(0,0,1):C.UVGC_721_262,(0,0,2):C.UVGC_721_263,(0,0,3):C.UVGC_721_264,(0,1,0):C.UVGC_720_259,(0,1,3):C.UVGC_720_260})

