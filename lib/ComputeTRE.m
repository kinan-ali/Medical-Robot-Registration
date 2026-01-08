function TRE = ComputeTRE()

global Tworld_base Tbase_eff Teff_inst Tworld_org torg_target

% Position actuelle du bout de l'instrument
Tworld_inst = Tworld_base * Tbase_eff * Teff_inst;
tworld_inst = Tworld_inst(1:3,4);

% Position de la cible
tworld_targ = Tworld_org(1:3,1:3) * torg_target + Tworld_org(1:3,4);

TRE = norm(tworld_targ - tworld_inst);