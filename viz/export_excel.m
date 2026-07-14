function export_excel(FIGS, xlsxpath)
%EXPORT_EXCEL Write one worksheet per figure, each containing the figure
%metadata (name, title, caption, explanation, interpretation) followed by
%every numeric-data block plotted in that figure.
%
% Purpose: Fulfils the requirement that the numeric data behind each
%          figure be exported to a single Excel workbook, one sheet per
%          figure, annotated with its caption, explanation and result
%          interpretation.
% Inputs : FIGS - struct array from make_figures.m
%          xlsxpath - full path of the .xlsx workbook to (re)create.
% Outputs: none (writes the workbook to disk).
% Reference: project post-processing requirements (CLAUDE.md).

if exist(xlsxpath,'file'), delete(xlsxpath); end

for i = 1:numel(FIGS)
    F = FIGS(i);
    sheet = sprintf('Figure_%03d', i);

    % --- header block (metadata) ---
    meta = { 'Figure name',        F.name; ...
             'Figure title',       F.title; ...
             'Caption',            char(F.caption); ...
             'Explanation',        char(F.explanation); ...
             'Interpretation',     char(F.interpretation) };
    writecell(meta, xlsxpath, 'Sheet', sheet, 'Range', 'A1');

    row = size(meta,1) + 2;     % leave a blank row

    % --- data blocks ---
    for b = 1:numel(F.blocks)
        B = F.blocks{b};
        writecell({['DATA BLOCK: ' B.label]}, xlsxpath, 'Sheet', sheet, 'Range', sprintf('A%d',row));
        row = row + 1;
        T = B.T;
        % include row names as a leading column when present
        if ~isempty(T.Properties.RowNames)
            rn = T.Properties.RowNames;
            T = addvars(T, rn, 'Before', 1, 'NewVariableNames', 'Name');
            T.Properties.RowNames = {};
        end
        writetable(T, xlsxpath, 'Sheet', sheet, 'Range', sprintf('A%d',row), 'WriteVariableNames', true);
        row = row + height(T) + 3;   % gap before next block
    end
end
fprintf('Excel export complete: %s (%d figure sheets)\n', xlsxpath, numel(FIGS));
end
