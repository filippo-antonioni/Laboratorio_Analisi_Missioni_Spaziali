clear
close all
clc

% =========================================================================
% OTTIMIZZATORE MULTI-OBIETTIVO (dV e TOF Ottimizzato)
% =========================================================================

% --- 1. CARICAMENTO DATI ---
T = load("DatiSC1-2026.txt");
mu = 398600;

gruppo = 23;

% Orbita finale/target
rr_f = T(gruppo,8:10)';
vv_f = T(gruppo,11:13)';
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% Orbita iniziale
orbita_iv = T(gruppo,2:7);
a_i  = orbita_iv(1);
e_i  = orbita_iv(2);
i_i  = orbita_iv(3);
OM_i = orbita_iv(4);
om_i = orbita_iv(5);
th_i = orbita_iv(6);

% --- 2. SETUP DELL'ANALISI OTTIMIZZAZIONE ---
tipi_bitangente = {'pa', 'ap', 'pp', 'aa'};
nomi_sequenze = {'A: Bitang -> Piano -> Peric', ...
                 'B: Piano -> Bitang -> Peric', ...
                 'C: Piano -> Peric -> Bitang (STANDARD)'};

risultati_Sequenza = strings(12,1);
risultati_Bitangente = strings(12,1);
risultati_dV_Bitangente = zeros(12,1);
risultati_dV_Piano = zeros(12,1);
risultati_dV_Pericentro = zeros(12,1);
risultati_dV_Totale = zeros(12,1);
risultati_TOF_Secondi = zeros(12,1);
risultati_TOF_Giorni = zeros(12,1); 

contatore = 1;

fprintf('--- INIZIO OTTIMIZZAZIONE MULTI-OBIETTIVO (dV e TOF) ---\n');
fprintf('Calcolo in corso per 12 combinazioni con ricerca del percorso minimo (Shortest Path)...\n\n');

% --- 3. CICLO DI OTTIMIZZAZIONE ---
for seq = 1:3
    for bit = 1:4
        tipo_b = tipi_bitangente{bit};
        
        % Determino anomalia vera di partenza e arrivo per il bitangente
        if strcmp(tipo_b, 'pa') || strcmp(tipo_b, 'pp')
            th_start_b = 0;
        else
            th_start_b = pi;
        end
        
        if strcmp(tipo_b, 'pa') || strcmp(tipo_b, 'aa')
            th_end_b = pi;
        else
            th_end_b = 0;
        end
        
        switch seq
            case 1 % A: Bitang -> Piano -> Peric
                [dV1_b, dV2_b, dt_b] = bitangentTransfer(a_i, e_i, a_f, e_f, tipo_b, mu);
                dV_b_tot = abs(dV1_b) + abs(dV2_b);
                
                [dV_pl, om_mid, th_node_base] = changeOrbitalPlane(a_f, e_f, i_i, OM_i, om_i, i_f, OM_f, mu);
                
                % [FIX TOF] Scelta del Nodo più vicino
                th_n1 = mod(th_node_base, 2*pi); th_n2 = mod(th_node_base + pi, 2*pi);
                dt1_1 = TOF(a_f, e_f, th_end_b, th_n1, mu);
                dt1_2 = TOF(a_f, e_f, th_end_b, th_n2, mu);
                if dt1_1 < dt1_2, dt1 = dt1_1; th_node = th_n1; else, dt1 = dt1_2; th_node = th_n2; end
                
                [dV_pe, thi_pe, thf_pe] = changePericenterArg(a_f, e_f, om_mid, om_f, mu);
                
                % [FIX TOF] Scelta dell'Anomalia più vicina
                dt2_1 = TOF(a_f, e_f, th_node, thi_pe(1), mu);
                dt2_2 = TOF(a_f, e_f, th_node, thi_pe(2), mu);
                if dt2_1 < dt2_2, dt2 = dt2_1; th_pe_f = thf_pe(1); else, dt2 = dt2_2; th_pe_f = thf_pe(2); end
                
                dt0 = TOF(a_i, e_i, th_i, th_start_b, mu);           
                dt3 = TOF(a_f, e_f, th_pe_f, th_f, mu);            
                TOF_tot = dt0 + dt_b + dt1 + dt2 + dt3;
                
            case 2 % B: Piano -> Bitang -> Peric
                [dV_pl, om_mid, th_node_base] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
                
                % [FIX TOF] Scelta Nodo
                th_n1 = mod(th_node_base, 2*pi); th_n2 = mod(th_node_base + pi, 2*pi);
                dt0_1 = TOF(a_i, e_i, th_i, th_n1, mu);
                dt0_2 = TOF(a_i, e_i, th_i, th_n2, mu);
                if dt0_1 < dt0_2, dt0 = dt0_1; th_node = th_n1; else, dt0 = dt0_2; th_node = th_n2; end
                
                [dV1_b, dV2_b, dt_b] = bitangentTransfer(a_i, e_i, a_f, e_f, tipo_b, mu);
                dV_b_tot = abs(dV1_b) + abs(dV2_b);
                
                [dV_pe, thi_pe, thf_pe] = changePericenterArg(a_f, e_f, om_mid, om_f, mu);
                
                % [FIX TOF] Scelta Pericentro
                dt2_1 = TOF(a_f, e_f, th_end_b, thi_pe(1), mu);
                dt2_2 = TOF(a_f, e_f, th_end_b, thi_pe(2), mu);
                if dt2_1 < dt2_2, dt2 = dt2_1; th_pe_f = thf_pe(1); else, dt2 = dt2_2; th_pe_f = thf_pe(2); end
                
                dt1 = TOF(a_i, e_i, th_node, th_start_b, mu);        
                dt3 = TOF(a_f, e_f, th_pe_f, th_f, mu);            
                TOF_tot = dt0 + dt1 + dt_b + dt2 + dt3;
                
            case 3 % C: Piano -> Pericentro -> Bitangente (STRATEGIA STANDARD)
                [dV_pl, om_mid, th_node_base] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
                
                % [FIX TOF] Scelta Nodo
                th_n1 = mod(th_node_base, 2*pi); th_n2 = mod(th_node_base + pi, 2*pi);
                dt0_1 = TOF(a_i, e_i, th_i, th_n1, mu);
                dt0_2 = TOF(a_i, e_i, th_i, th_n2, mu);
                if dt0_1 < dt0_2, dt0 = dt0_1; th_node = th_n1; else, dt0 = dt0_2; th_node = th_n2; end
                
                [dV_pe, thi_pe, thf_pe] = changePericenterArg(a_i, e_i, om_mid, om_f, mu);
                
                % [FIX TOF] Scelta Pericentro
                dt1_1 = TOF(a_i, e_i, th_node, thi_pe(1), mu);
                dt1_2 = TOF(a_i, e_i, th_node, thi_pe(2), mu);
                if dt1_1 < dt1_2, dt1 = dt1_1; th_pe_f = thf_pe(1); else, dt1 = dt1_2; th_pe_f = thf_pe(2); end
                
                [dV1_b, dV2_b, dt_b] = bitangentTransfer(a_i, e_i, a_f, e_f, tipo_b, mu);
                dV_b_tot = abs(dV1_b) + abs(dV2_b);
                
                dt2 = TOF(a_i, e_i, th_pe_f, th_start_b, mu);      
                dt3 = TOF(a_f, e_f, th_end_b, th_f, mu);             
                TOF_tot = dt0 + dt1 + dt2 + dt_b + dt3;
        end
        
        dV_Totale = dV_b_tot + abs(dV_pl) + abs(dV_pe);
        
        risultati_Sequenza(contatore) = nomi_sequenze{seq};
        risultati_Bitangente(contatore) = string(tipo_b);
        risultati_dV_Bitangente(contatore) = dV_b_tot;
        risultati_dV_Piano(contatore) = abs(dV_pl);
        risultati_dV_Pericentro(contatore) = abs(dV_pe);
        risultati_dV_Totale(contatore) = dV_Totale;
        risultati_TOF_Secondi(contatore) = TOF_tot;
        risultati_TOF_Giorni(contatore) = TOF_tot / (24*3600);
        
        contatore = contatore + 1;
    end
end

% --- 4. CREAZIONE E ORDINAMENTO DELLA TABELLA ---
TabellaRisultati = table(risultati_Sequenza, risultati_Bitangente, ...
    risultati_dV_Bitangente, risultati_dV_Piano, risultati_dV_Pericentro, risultati_dV_Totale, risultati_TOF_Giorni, ...
    'VariableNames', {'Sequenza', 'Tipo_Bitang', 'dV_Bit', 'dV_Piano', 'dV_Peric', 'dV_TOT', 'TOF_Giorni'});

TabellaOrdinata_dV = sortrows(TabellaRisultati, 'dV_TOT');
disp('CLASSIFICA COMPLETA ORDINATA PER DELTA V:');
disp(TabellaOrdinata_dV);

migliore = TabellaOrdinata_dV(1, :);
fprintf('\n*******************************************************\n');
fprintf('LA STRATEGIA OTTIMA ASSOLUTA E'':\n');
fprintf('Sequenza:         %s\n', migliore.Sequenza);
fprintf('Tipo Bitangente:  %s\n', migliore.Tipo_Bitang);
fprintf('Costo Totale:     %.4f km/s\n', migliore.dV_TOT);
fprintf('Tempo Richiesto:  %.2f Giorni\n', migliore.TOF_Giorni);
fprintf('*******************************************************\n\n');


% =========================================================================
% 5. PLOT 3D: STRATEGIA STANDARD (Piano -> Pericentro -> Bitangente 'pa')
% =========================================================================
% Ricalcolo forzato della cinematica per la sequenza Standard
[~, om_mid, th_node_base] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);

th_n1 = mod(th_node_base, 2*pi); 
th_n2 = mod(th_node_base + pi, 2*pi);
dt0_1 = TOF(a_i, e_i, th_i, th_n1, mu);
dt0_2 = TOF(a_i, e_i, th_i, th_n2, mu);
if dt0_1 < dt0_2, th_node = th_n1; else, th_node = th_n2; end

[~, thi_pe, thf_pe] = changePericenterArg(a_i, e_i, om_mid, om_f, mu);

dt1_1 = TOF(a_i, e_i, th_node, thi_pe(1), mu);
dt1_2 = TOF(a_i, e_i, th_node, thi_pe(2), mu);
if dt1_1 < dt1_2, th_pe_start = thi_pe(1); else, th_pe_start = thi_pe(2); end

% Per il bitangente 'pa', partenza dal pericentro (0) e arrivo all'apocentro (pi)
th_start_b = 0; 
th_end_b = pi;

r1_mag = a_i*(1-e_i);
r2_mag = a_f*(1+e_f);
a_t = (r1_mag + r2_mag) / 2;
e_t = abs(r2_mag - r1_mag) / (r1_mag + r2_mag);

% --- SETUP FIGURA ---
fig_std = figure('Name', 'Strategia Standard (Piano -> Pericentro -> Bitangente pa)', 'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(45, 30); rotate3d on;
xlabel('X [km]', 'Interpreter', 'latex', 'FontSize', 12); 
ylabel('Y [km]', 'Interpreter', 'latex', 'FontSize', 12); 
zlabel('Z [km]', 'Interpreter', 'latex', 'FontSize', 12);
title('\textbf{Strategia Standard: Piano $\rightarrow$ Pericentro $\rightarrow$ Bitangente ''pa''}', 'Interpreter', 'latex', 'Color', 'w', 'FontSize', 15);

% Terra Fotorealistica
[xE, yE, zE] = sphere(50);
try
    load topo topo topomap1;
    surf(xE * 6371, yE * 6371, zE * 6371, 'FaceColor', 'texturemap', 'CData', topo, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(ax, topomap1);
catch
    surf(xE * 6371, yE * 6371, zE * 6371, 'FaceColor', [0.1 0.4 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

th_vec_plot = linspace(0, 2*pi, 300);

% --- PLOT ORBITE (Colori Immagine) ---
% 1. Iniziale (Blu)
plotOrbit_Fast(a_i, e_i, i_i, OM_i, om_i, mu, th_vec_plot, [0 0.447 0.741], '-', 1.5, 'Orbita Iniziale');

% 2. Post-Cambio Piano (Magenta)
plotOrbit_Fast(a_i, e_i, i_f, OM_f, om_mid, mu, th_vec_plot, [0.494 0.184 0.556], '-', 1.5, 'Parcheggio (Post-Piano)');

% 3. Post-Cambio Pericentro (Ciano)
plotOrbit_Fast(a_i, e_i, i_f, OM_f, om_f, mu, th_vec_plot, [0.301 0.745 0.933], '-', 1.5, 'Parcheggio (Post-Pericentro)');

% 4. Trasferimento Bitangente 'pa' (Rosso Tratteggiato Spesso)
th_vec_trasf = linspace(th_start_b, th_end_b, 150);
plotOrbit_Fast(a_t, e_t, i_f, OM_f, om_f, mu, th_vec_trasf, [0.635 0.078 0.184], '--', 2.5, 'Trasferimento Bitangente (pa)');

% 5. Finale Target (Verde Spessa)
plotOrbit_Fast(a_f, e_f, i_f, OM_f, om_f, mu, th_vec_plot, [0.466 0.674 0.188], '-', 2.5, 'Orbita Finale Target');

% --- MARKERS (Esatti da Immagine) ---
% M1: Cambio Piano (Triangolo Magenta)
[r_m1,~] = par2car(a_i, e_i, i_i, OM_i, om_i, th_node, mu);
plot3(r_m1(1), r_m1(2), r_m1(3), '^k', 'MarkerSize', 11, 'MarkerFaceColor', [0.494 0.184 0.556], 'DisplayName', 'M1: Cambio Piano');

% M2: Cambio Pericentro (Stella Ciano)
[r_m2,~] = par2car(a_i, e_i, i_f, OM_f, om_mid, th_pe_start, mu);
plot3(r_m2(1), r_m2(2), r_m2(3), 'pk', 'MarkerSize', 11, 'MarkerFaceColor', [0.301 0.745 0.933], 'DisplayName', 'M2: Cambio Pericentro');

% M3: Inizio Bitangente (Rombo Rosso)
[r_m3,~] = par2car(a_i, e_i, i_f, OM_f, om_f, th_start_b, mu);
plot3(r_m3(1), r_m3(2), r_m3(3), 'dk', 'MarkerSize', 9, 'MarkerFaceColor', [0.635 0.078 0.184], 'DisplayName', 'M3: Inizio Bitangente');

% M4: Arrivo Target (Stella Verde)
[r_m4,~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_end_b, mu);
plot3(r_m4(1), r_m4(2), r_m4(3), 'pk', 'MarkerSize', 13, 'MarkerFaceColor', [0.466 0.674 0.188], 'MarkerEdgeColor', 'w', 'DisplayName', 'M4: Arrivo Target');

legend('show', 'Location', 'bestoutside', 'FontSize', 10, 'TextColor', 'w', 'Color', 'k', 'EdgeColor', [0.5 0.5 0.5], 'Interpreter', 'latex');
hold off;

% =========================================================================
% FUNZIONI AUSILIARIE LOCALI
% =========================================================================
function plotOrbit_Fast(a, e, i, OM, om, mu, th_vec, col, stile, width, nome)
    pts = zeros(3, length(th_vec));
    for j = 1:length(th_vec)
        [pts(:,j), ~] = par2car(a, e, i, OM, om, th_vec(j), mu);
    end
    plot3(pts(1,:), pts(2,:), pts(3,:), 'Color', col, 'LineStyle', stile, ...
          'LineWidth', width, 'DisplayName', ['\textbf{', nome, '}']);
end
