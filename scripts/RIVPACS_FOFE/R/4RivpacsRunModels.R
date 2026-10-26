# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #### Hessun oma final mallin ajo 23.5.2023  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aja Sicklen RF malli
# source("D:/Aroviita/3Rproj/1Functions/predictive_model_scripts_v4.2.1/model.predict.RanFor.4.2.r")
# source("M:/GVESIBIO/3Rproj/1Functions/predictive_model_scripts_v4.2.1/model.predict.RanFor.4.2.r")
source("R/model.predict.RanFor.4.2.r")

# test
# OERFtest <- model.predict.RanFor.4.2(bugcal.pa = BIO.cal.pa, grps.final = grps.finalAB2, preds.final = preds.finalAB2, ranfor.mod = rf.modAB, prednew = ENVtest, bugnew = BIOtest, Pc = 0.5, Cal.OOB = FALSE)
# # Tee kuva
# pdf(paste(sdirfigs, paste("RFOEabuscatter_N", nrow(OERFtest$OE.scores), ".pdf", sep = ""), sep = ""))
# plot(OERFtest$OE.scores[c("O", "E", "OoverE", "BC")], main = paste("Final RF-model of ", nrow(OERFtest$OE.scores), " TEST sites", sep = ""))
# dev.off()

# vld
# OERFvld <- model.predict.RanFor.4.2(bugcal.pa = BIO.cal.pa, grps.final = grps.finalAB2, preds.final = preds.finalAB2, ranfor.mod = rf.modAB, prednew = ENVvld, bugnew = BIOvld, Pc=0.5,Cal.OOB=FALSE)
# # Tee kuva
# pdf(paste(sdirfigs, paste("RFOEabuscatter_N", nrow(OERFvld$OE.scores), ".pdf", sep=""), sep=""))
# plot(OERFvld$OE.scores[c("O","E","OoverE", "BC")], main=paste("Final RF-model of ", nrow(OERFvld$OE.scores), " VLD sites", sep=""))
# dev.off()

# all
OERFall <- model.predict.RanFor.4.2(bugcal.pa = BIO.cal.pa, grps.final = grps.finalAB2, preds.final = preds.finalAB2, ranfor.mod = rf.modAB, prednew = ENVall, bugnew = BIOall, Pc = 0.5, Cal.OOB = FALSE)
# Tee kuva
pdf(paste(sdirfigs, paste("RFOEabuscatter_N", nrow(OERFall$OE.scores), ".pdf", sep = ""), sep = ""))
plot(OERFall$OE.scores[c("O", "E", "OoverE", "BC")], main = paste("Final RF-model of ", nrow(OERFall$OE.scores), " ALL sites", sep = ""))
dev.off()


# Tallenna
# write.table(OERFvld$OE.scores, sep=";", col.names=NA, file=paste(sdir, "RF_OE.csv", sep=""))
# write.table(OERFvld$Capture.Probs, sep=";", col.names=NA, file=paste(sdir, "RF_CaptP.csv", sep=""))
# write.table(OERFvld$Group.Occurrence.Probs, sep=";", col.names=NA, file=paste(sdir, "RF_GrpOP.csv", sep=""))


# Lis?? OE:t envidataan ja tallenna
# On sama j?rjestys
# po_env_N648$FH_id_uniikki == row.names(OERF$OE.scores)
# (po_env_N648$OJITETUR, d$RFOE)
# po_env_N648 <- cbind(po_env_N648, OERF$OE.scores)
# fwrite(po_env_N648, sep=";", paste(sdir, "po_env_N648.csv", sep=""))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# T?m? skripti Freshabitin purodatan:
# Abundanssimalli RF
# source("D:/Aroviita/3Rproj/1Functions/RIVPACS_Abundance/model.predict.RanFor.4.2_and_abuOE.r");
# marian malli
# OE.assess.IMP2_abu_<-model.predict.RanFor.4.2_abu(bugcal.abu,bugcal.pa,grps.final_5,preds.final_REF,ranfor.mod=rf.mod_kal_5, prednew=pred_IMP2,bugnew=bug_IMP2,Pc=0.5,Cal.OOB=FALSE);

# koko N=648 -datalla
# OEaRF <- model.predict.RanFor.4.2_abu(bugcal.abu=bug_pohjis_REF, bugcal.pa=bug_pohjis_REF.pa, grps.final=grps.final_5, preds.final=preds.final_REF, ranfor.mod=rf.mod_kal_5, prednew=pred_po, bugnew=po_spedf, Pc=0.5, Cal.OOB=FALSE);

# Tee kuva
# pdf(paste(sdirfigs, "RFOEvsBCvsBrayabu_scatter.pdf", sep=""))
# plot(OEaRF$OE.scores[c("OoverE", "BC", "Bray_OEabu", "Bray_OEabulog")], main="Final RF-model abundance invertebrates all 648 sites")
# dev.off()


# Tallenna
# write.table(OEaRF$OE.scores, sep=";", col.names=NA, file=paste(sdir, "RFabu_OE_N648.csv", sep=""))
# write.table(OEaRF$Capture.Probs, sep=";", col.names=NA, file=paste(sdir, "RFabu_CaptP_N648.csv", sep=""))
# write.table(OEaRF$Group.Occurrence.Probs, sep=";", col.names=NA, file=paste(sdir, "RFabu_GrpOP_N648.csv", sep=""))
# write.table(OEaRF$Exp.Abus, sep=";", col.names=NA, file=paste(sdir, "RFabu_ExpAbus_N648.csv", sep=""))
# write.table(OEaRF$Exp.LogAbus, sep=";", col.names=NA, file=paste(sdir, "RFabu_ExpLogAbus_N648.csv", sep=""))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Show all items
print(ls())
