#######Mean cluster generation 

#done
d[, "Glucose_g"] <- rowMeans(d[, c("glc", "glucose")])
#done
d[, "s-HDL_g"] <- rowMeans(d[, c("shdlp", "shdlpl", "shdlfc", "shdll")])
#done
d[, "Val-Leu_g"] <- rowMeans(d[, c("val", "leu")])

#done
d[, "m/xl-VLDL-C_g"] <- rowMeans(d[, c("xxlvldlce", "mvldlc", "mvldlce")])

#done
d[, "s/Tot-VLDL-C_g"] <- rowMeans(d[, c("apob", "remnantc", "vldlc", "svldlc", "svldlce" )])
#done
d[, "(x)(x)l-VLDL_g"] <- rowMeans(d[, c("vldld", "xxlvldlc", "xlvldlc", "xlvldlce", "xlvldlfc", "xlvldlpl", "xxlvldlp", "xxlvldll", "xxlvldltg", "xxlvldlfc", "xxlvldlpl", "serumtg", "svldltg", "lvldlc", "lvldlce", "mvldll", "mvldlp", "lvldlfc", "lvldlpl", "mvldlfc", "mvldlpl", "xlvldltg", "xlvldll", "xlvldlp", "mvldltg", "vldltg", "lvldltg", "lvldll", "lvldlp")])

#done
d[, "HDL-TG_g"] <- rowMeans(d[, c("shdltg", "hdltg", "mhdltg")])
#done
d[, "s-VLDL_g"] <- rowMeans(d[, c("mufa", "xsvldltg", "svldll", "svldlp", "svldlfc", "svldlpl")])

#done
d[, "xl-HDL-C_g"] <- rowMeans(d[, c("xlhdlce", "xlhdlc")])

#done
d[, "Apoa1_&_HDL3-C_g"] <- rowMeans(d[, c("apoa1", "hdl3c")])

#done
d[, "l-HDL_g"] <- rowMeans(d[, c("hdl2c", "hdlc", "lhdlpl", "lhdll", "lhdlp", "lhdlfc", "lhdlc", "lhdlce", "xlhdlfc", "xlhdll", "hdld", "xlhdlp", "xlhdlpl")])

#done
d[, "s-HDL-C_g"] <- rowMeans(d[, c("shdlc", "shdlce")])


d[, "m-HDL_g"] <- rowMeans(d[, c("mhdlpl", "mhdll", "mhdlp", "mhdlfc", "mhdlc", "mhdlce")])

#done
d[, "xs-VLDL-C_g"] <- rowMeans(d[, c("xsvldlc", "xsvldlce")])

#done
d[, "xs-VLDL_g"] <- rowMeans(d[, c("xsvldlfc", "xsvldll", "xsvldlp")])
#done
d[, "PUFA_g"] <- rowMeans(d[, c("la", "faw6", "pufa")])
#done
d[, "LDL_&_IDL_g"] <- rowMeans(d[, c("mldlp", "lldlp", "mldll", "sldll", "sldlp", "estc", "serumc", "mldlpl", "sldlpl", "freec", "ldlc", "lldlc", "lldlce", "lldll", "lldlpl", "mldlc", "mldlce", "sldlc", "sldlce", "idll", "idlp", "idlpl", "idlfc", "lldlfc", "mldlfc", "sldlfc", "xsvldlpl", "idlc", "idlce")])
#done
d[, "LDL-TG_g"] <- rowMeans(d[, c("sldltg", "idltg", "mldltg", "ldltg", "lldltg")])

d[, "l-HDL-TG_g"] <- rowMeans(d[, c("lhdltg", "xlhdltg")])
#done
d[, "SFA_&_Tot-FA_g"] <- rowMeans(d[, c("sfa", "totfa")])
#done
d[, "Choline_g"] <- rowMeans(d[, c("totpg", "totcho", "pc")])


