'use strict';

$.getJSON('/static/phenotypes.json').then(function(data) {
    $(document).ready(function() {
        var table = new Tabulator('#table', {
//            height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local',
            paginationSize: 100,
            columns: [
                {title: 'Code', field:'phecode', formatter:'link', formatterParams: { urlPrefix: '/pheno/' }, headerFilter:true},
                {title: 'Name', field:'phenostring', formatter:'link', formatterParams: {url: function(cell){return '/pheno/'+cell.getData().phecode;}}, headerFilter:true, widthGrow:2},
                {title: 'Category', field:'category', headerFilter:true, widthGrow:1},
                {title: 'Sample Size', field:'sample_size', formatter:'comma_fmt'},
                //{title: '#Controls', field:'num_controls', formatter:'comma_fmt'},
                //{title: 'Excluded phecodes', field:'phecode_exclude_range', formatter:function(cell){return cell.getValue()+' ('+cell.getData().phecode_exclude_description+')'}},
                //{title: 'Sex', field:'sex'},
                {title: 'Num p&le;2.5<sup>-6</sup> Associations', field:'num_sig_assocs', width:200},
                {title: 'Best Gene', field:'best_assoc'},
                {title: 'Best P-value', field:'best_pval', formatter:'2digit_fmt'},
            ],
            data: data,
            initialSort: [{column:"best_pval", dir:"asc"}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});
