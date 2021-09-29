PtrList<radiation::radiationModel> radiation1(2);
PtrList<radiation::radiationModel> radiation2(2);
PtrList<radiation::radiationModel> radiation3(2);

radiation1.set(1,radiation::radiationModel::New(mixture.T()));
radiation2.set(1,radiation::radiationModel::New(mixture.T()));
radiation3.set(1,radiation::radiationModel::New(mixture.T()));

radiation::radiationModel& rad1 = radiation1[1];
radiation::radiationModel& rad2 = radiation2[1];
radiation::radiationModel& rad3 = radiation3[1];