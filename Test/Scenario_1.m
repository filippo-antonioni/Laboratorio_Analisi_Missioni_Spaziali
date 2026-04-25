clear
close all
clc

% --- 1. CARICAMENTO DATI ---
T = load("DatiSC1-2026.txt");
mu = 398600;

% Orbita target/finale (vettori di stato)
rr_f = T(23,8:10)';
vv_f = T(23,11:13)';
state_f = [rr_f; vv_f];

% Trasformazione per avere i parametri kepleriani bersaglio
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% Orbita iniziale (parametri kepleriani dal file)
orbita_iv = T(23,2:7);
a_i  = orbita_iv(1);
e_i  = orbita_iv(2);
i_i  = orbita_iv(3);
OM_i = orbita_iv(4);
om_i = orbita_iv(5);
th_i = orbita_iv(6);

state_i = [a_i e_i i_i OM_i om_i th_i]';

% --- 2. CALCOLO MANOVRE ---
delta_v = [];
delta_t = [];
type_v = ['ap'; 'pa'; 'aa'; 'pp'];




% 2.1 Calcolo i costi di tutti i trasferimenti bitangenti
for k = 1:4 
    type = type_v(k,:);
    [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu);
    % Sommo i valori assoluti perché rappresentano il propellente consumato
    delta_v = [delta_v; abs(DeltaV1) + abs(DeltaV2)]; 
    delta_t = [delta_t; Deltat];
end



% SELEZIONA IL TIPO DI MANOVRA BITANGENTE DA VISUALIZZARE
% (1='ap', 2='pa', 3='aa', 4='pp')
indice_scelto = 2; 
tipo_scelto = type_v(indice_scelto, :);

deltav_bitang = delta_v(indice_scelto);
deltat_bitang = delta_t(indice_scelto);

% 2.2 Cambio di Piano 
% Ottiene om_fun, la nuova anomalia del pericentro forzata dalla geometria sferica
[deltav_plane, om_fun, theta] = changeOrbitalPlane(a_f, e_f, i_i, OM_i, om_i, i_f, OM_f, mu);

% 2.3 Cambio Anomalia Pericentro 
% Ruoto il pericentro partendo da om_fun fino al bersaglio finale om_f
[deltav_arg, thi_fun, thf_fun] = changePericenterArg(a_f, e_f, om_fun, om_f, mu);

% Costo totale della sequenza
DeltaV_Tot = deltav_bitang + deltav_plane + deltav_arg;

fprintf('\n--- RISULTATI MANOVRA (Sequenza Bitangente: %s) ---\n', tipo_scelto);
fprintf('DeltaV Bitangente:        %.4f km/s\n', deltav_bitang);
fprintf('DeltaV Cambio Piano:      %.4f km/s\n', deltav_plane);
fprintf('DeltaV Cambio Pericen.:   %.4f km/s\n', deltav_arg);
fprintf('DELTAV TOTALE:            %.4f km/s\n\n', DeltaV_Tot);


% --- 4. CALCOLO DEL TEMPO DI VOLO TOTALE (Fasi di attesa) ---

% 1. Determiniamo l'anomalia vera in cui arriviamo dopo il bitangente
% (Se partiamo per l'apocentro arriviamo a pi, se partiamo per il pericentro arriviamo a 0)
switch tipo_scelto
    case 'pa'
        th_arrivo_bitang=pi;
        th_inizio_bitang=0;
    case 'aa'
        th_arrivo_bitang=pi;
        th_inizio_bitang=pi;
    case 'ap'
        th_arrivo_bitang = 0;
        th_inizio_bitang=pi;
    case 'pp'
        th_arrivo_bitang = 0;
        th_inizio_bitang=0;
    otherwise
        error('Errore nel type della manovra bitangente');

end


% if strcmp(tipo_scelto, 'pa') || strcmp(tipo_scelto, 'aa')
%     th_arrivo_bitang = pi;
% else
%     th_arrivo_bitang = 0;
% end

% calcolo tempo coasting da th_i a punto in cui faccio il bitang
deltat_coast1=TOF(a_i,e_i,th_i,th_inizio_bitang,mu);

% 2. Fase di attesa (volo passivo) dall'arrivo del bitangente fino al nodo del cambio piano
deltat_attesa1 = TOF(a_f, e_f, th_arrivo_bitang, theta, mu);

% 3. Fase di attesa dal cambio piano fino alla rotazione del pericentro
% Usiamo la prima intersezione utile fornita dalla funzione (thi_fun(1))
deltat_attesa2 = TOF(a_f, e_f, theta, thi_fun(1), mu);

% tempo per arrivare da theta cambio pericentro fino a theta finale
deltat_coast3=TOF(a_f,e_f,thf_fun(1),th_f,mu);

% TEMPO DI VOLO TOTALE
TOF_Totale = deltat_coast1+deltat_bitang + deltat_attesa1 + deltat_attesa2 +deltat_coast3;
TOF_Totale_giorni = TOF_Totale / (24*3600); % Conversione in giorni

% Aggiungiamo le stampe a schermo
fprintf('--- RISULTATI TEMPI (Volo Passivo e Attivo) ---\n');
fprintf('Trasferimento Bitangente: %.2f s (%.2f giorni)\n', deltat_bitang, deltat_bitang/(24*3600));
fprintf('Attesa 1 (verso il nodo): %.2f s (%.2f giorni)\n', deltat_attesa1, deltat_attesa1/(24*3600));
fprintf('Attesa 2 (verso manovra): %.2f s (%.2f giorni)\n', deltat_attesa2, deltat_attesa2/(24*3600));
fprintf('TEMPO DI VOLO TOTALE:     %.2f s (%.2f giorni)\n\n', TOF_Totale, TOF_Totale_giorni);


% --- 3. GRAFICA AVANZATA ---
% Generiamo la scena base con plotTransfer 
% (Disegna: Terra, Orbita Iniziale blu, Orbita Finale verde e archi di manovra)
plotTransfer(state_i, 'par', state_f, 'car', tipo_scelto, mu);

hold on;
dth = 0.05;
th_vec = 0:dth:2*pi;

% --- AGGIUNTA 1: Orbita post-Bitangente (Linea continua) ---
% Ha la forma e dimensione finale (a_f, e_f) ma il piano (i_i, OM_i) 
% e il pericentro (om_i) sono ancora quelli dell'orbita iniziale.
rr_post_bitang = zeros(3, length(th_vec));
for k = 1:length(th_vec)
    [rr_k, ~] = par2car(a_f, e_f, i_i, OM_i, om_i, th_vec(k), mu);
    rr_post_bitang(:, k) = rr_k;
end
plot3(rr_post_bitang(1,:), rr_post_bitang(2,:), rr_post_bitang(3,:), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Orbita post-Bitangente');

% --- AGGIUNTA 2: Orbita post-Cambio Piano (Linea tratteggiata) ---
% Ha forma a_f, e_f e piano finale i_f, OM_f, ma il pericentro è momentaneamente om_fun.
rr_intermedia = zeros(3, length(th_vec));
for k = 1:length(th_vec)
    [rr_k, ~] = par2car(a_f, e_f, i_f, OM_f, om_fun, th_vec(k), mu);
    rr_intermedia(:, k) = rr_k;
end
plot3(rr_intermedia(1,:), rr_intermedia(2,:), rr_intermedia(3,:), 'c-.', 'LineWidth', 1.5, 'DisplayName', 'Orbita post-Cambio Piano');

% --- AGGIUNTA 3: MARKER DEI PUNTI DI MANOVRA ---
% Utilizziamo par2car per trovare le coordinate 3D dei punti in cui
% il satellite accende i motori, usando i parametri orbitali corretti per quella fase.

% 1. Inizio Bitangente (Avviene sull'orbita iniziale, ad anomalia th_inizio_bitang)
[r_man1, ~] = par2car(a_i, e_i, i_i, OM_i, om_i, th_inizio_bitang, mu);
plot3(r_man1(1), r_man1(2), r_man1(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', '1: Inizio Bitangente');

% 2. Fine Bitangente (Avviene sull'orbita post-bitangente, ad anomalia th_arrivo_bitang)
[r_man2, ~] = par2car(a_f, e_f, i_i, OM_i, om_i, th_arrivo_bitang, mu);
plot3(r_man2(1), r_man2(2), r_man2(3), 's', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'DisplayName', '2: Fine Bitangente');

% 3. Cambio di Piano (Avviene sull'orbita post-bitangente, ad anomalia theta)
[r_man3, ~] = par2car(a_f, e_f, i_i, OM_i, om_i, theta, mu);
plot3(r_man3(1), r_man3(2), r_man3(3), '^', 'MarkerSize', 9, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'DisplayName', '3: Cambio Piano');

% 4. Cambio Pericentro (Avviene sull'orbita post-cambio piano, ad anomalia thi_fun(1))
[r_man4, ~] = par2car(a_f, e_f, i_f, OM_f, om_fun, thi_fun(1), mu);
plot3(r_man4(1), r_man4(2), r_man4(3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', '4: Cambio Pericentro');
% Aggiorno la visualizzazione
legend('show', 'Location', 'best');
view(3);
hold off;

