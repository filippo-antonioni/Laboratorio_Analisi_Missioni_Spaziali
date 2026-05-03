clear
close all
clc

%  Ottimizzatore 

% --- 1. CARICAMENTO DATI ---
T = load("DatiSC1-2026.txt");
mu = 398600;

gruppo=23;

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
                 'C: Piano -> Peric -> Bitang'};

% Pre-alloco variabili per la tabella (usando il comando strings raccomandato)
risultati_Sequenza = strings(12,1);
risultati_Bitangente = strings(12,1);
risultati_dV_Bitangente = zeros(12,1);
risultati_dV_Piano = zeros(12,1);
risultati_dV_Pericentro = zeros(12,1);
risultati_dV_Totale = zeros(12,1);
risultati_TOF_Secondi = zeros(12,1);
risultati_TOF_Giorni = zeros(12,1); % Aggiungo in giorni per leggibilità

contatore = 1;

fprintf('--- INIZIO OTTIMIZZAZIONE MULTI-OBIETTIVO (dV e TOF) ---\n');
fprintf('Calcolo in corso per 12 combinazioni...\n\n');

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
            case 1
                % =========================================================
                % SEQUENZA A: Bitangente -> Piano -> Pericentro
                % =========================================================
                [dV1_b, dV2_b, dt_b] = bitangentTransfer(a_i, e_i, a_f, e_f, tipo_b, mu);
                dV_b_tot = abs(dV1_b) + abs(dV2_b);
                
                [dV_pl, om_mid, th_node] = changeOrbitalPlane(a_f, e_f, i_i, OM_i, om_i, i_f, OM_f, mu);
                
                [dV_pe, thi_pe, thf_pe] = changePericenterArg(a_f, e_f, om_mid, om_f, mu);
                
                % Calcolo Tempi Sequenza A
                dt0 = TOF(a_i, e_i, th_i, th_start_b, mu);           % Attesa pre-bitangente
                dt1 = TOF(a_f, e_f, th_end_b, th_node, mu);          % Attesa post-bitangente fino a nodo
                dt2 = TOF(a_f, e_f, th_node, thi_pe(1), mu);         % Attesa da nodo a rotazione
                dt3 = TOF(a_f, e_f, thf_pe(1), th_f, mu);            % Attesa da rotazione a target
                TOF_tot = dt0 + dt_b + dt1 + dt2 + dt3;
                
            case 2
                % =========================================================
                % SEQUENZA B: Piano -> Bitangente -> Pericentro
                % =========================================================
                [dV_pl, om_mid, th_node] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
                
                [dV1_b, dV2_b, dt_b] = bitangentTransfer(a_i, e_i, a_f, e_f, tipo_b, mu);
                dV_b_tot = abs(dV1_b) + abs(dV2_b);
                
                [dV_pe, thi_pe, thf_pe] = changePericenterArg(a_f, e_f, om_mid, om_f, mu);
                
                % Calcolo Tempi Sequenza B
                dt0 = TOF(a_i, e_i, th_i, th_node, mu);              % Attesa fino a nodo cambio piano
                dt1 = TOF(a_i, e_i, th_node, th_start_b, mu);        % Attesa da nodo a bitangente
                dt2 = TOF(a_f, e_f, th_end_b, thi_pe(1), mu);        % Attesa post-bitangente fino a rotazione
                dt3 = TOF(a_f, e_f, thf_pe(1), th_f, mu);            % Attesa da rotazione a target
                TOF_tot = dt0 + dt1 + dt_b + dt2 + dt3;
                
            case 3
                % =========================================================
                % SEQUENZA C: Piano -> Pericentro -> Bitangente
                % =========================================================
                [dV_pl, om_mid, th_node] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);
                
                [dV_pe, thi_pe, thf_pe] = changePericenterArg(a_i, e_i, om_mid, om_f, mu);
                
                [dV1_b, dV2_b, dt_b] = bitangentTransfer(a_i, e_i, a_f, e_f, tipo_b, mu);
                dV_b_tot = abs(dV1_b) + abs(dV2_b);
                
                % Calcolo Tempi Sequenza C
                dt0 = TOF(a_i, e_i, th_i, th_node, mu);              % Attesa fino a nodo
                dt1 = TOF(a_i, e_i, th_node, thi_pe(1), mu);         % Attesa da nodo a rotazione pericentro
                dt2 = TOF(a_i, e_i, thf_pe(1), th_start_b, mu);      % Attesa da rotazione a bitangente
                dt3 = TOF(a_f, e_f, th_end_b, th_f, mu);             % Attesa post-bitangente fino a target
                TOF_tot = dt0 + dt1 + dt2 + dt_b + dt3;
        end
        
        % Calcolo del totale dV
        dV_Totale = dV_b_tot + abs(dV_pl) + abs(dV_pe);
        
        % Salvataggio nei vettori
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
% Inserisco tutto in una tabella (ho abbreviato i nomi delle variabili dV per compattezza di visualizzazione)
TabellaRisultati = table(risultati_Sequenza, risultati_Bitangente, ...
    risultati_dV_Bitangente, risultati_dV_Piano, risultati_dV_Pericentro, risultati_dV_Totale, risultati_TOF_Giorni, ...
    'VariableNames', {'Sequenza', 'Tipo_Bitang', 'dV_Bit', 'dV_Piano', 'dV_Peric', 'dV_TOT', 'TOF_Giorni'});

% Riordino la tabella in ordine crescente basandomi PRIMA sul Delta V Totale
TabellaOrdinata_dV = sortrows(TabellaRisultati, 'dV_TOT');

% --- 5. STAMPA A SCHERMO DEI RISULTATI ---
disp('CLASSIFICA COMPLETA ORDINATA PER DELTA V (Propellente Minimo):');
disp(TabellaOrdinata_dV);

% Estrazione della migliore per Delta V
migliore = TabellaOrdinata_dV(1, :);

fprintf('\n*******************************************************\n');
fprintf('LA STRATEGIA OTTIMA ASSOLUTA (Minimo Propellente) E'':\n');
fprintf('Sequenza:         %s\n', migliore.Sequenza);
fprintf('Tipo Bitangente:  %s\n', migliore.Tipo_Bitang);
fprintf('Costo Totale:     %.4f km/s\n', migliore.dV_TOT);
fprintf('Tempo Richiesto:  %.2f Giorni\n', migliore.TOF_Giorni);
fprintf('*******************************************************\n\n');

% Bonus: Trova la più veloce in assoluto (magari costa troppo in Delta V, ma è interessante da vedere)
TabellaOrdinata_TOF = sortrows(TabellaRisultati, 'TOF_Giorni');
piu_veloce = TabellaOrdinata_TOF(1, :);

fprintf('Nota sulla rapidità:\n');
fprintf('La missione piu'' veloce in assoluto dura %.2f giorni (Sequenza %s con bitangente %s),\n', piu_veloce.TOF_Giorni, piu_veloce.Sequenza, piu_veloce.Tipo_Bitang);
fprintf('ma richiede un Delta V di %.4f km/s.\n', piu_veloce.dV_TOT);