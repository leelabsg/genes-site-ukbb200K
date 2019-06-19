'use strict';

$.getJSON('/api/variants/'+model.genename+'/'+model.phecode).then(function(resp) {
    window._debug.resp = resp;
    var df = resp.df;
    df.pos = df.pos_delta;
    delete df.pos_delta;
    _.map(df.pos, function(elem, i) { if (i>0) { df.pos[i] += df.pos[i-1]; } })
    var data = dataframe_to_objects(df);
    console.log(data);

    $(function() {
        var table = new Tabulator('#table', {
            //height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local', // TODO: disable pagination if <100 variants
            paginationSize: 100,
            columns: [
                {title: 'Position on chr'+resp.chrom, field:'pos'},
                {title: 'Allele', field:'base'},
                {title: 'MAF (Minor Allele Frequency)', field:'maf'},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
                {title: 'P-value', field:'pval'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipsHeader:true,
            tooltipGenerationMode: 'hover', // generate tooltips just-in-time when the data is hovered
            tooltips: function(cell) {
                // this function attempts to check whether an ellipsis ('...') is hiding part of the data.
                // to do so, I compare element.clientWidth against element.scrollWidth;
                // when scrollWidth is bigger, that means we're hiding part of the data.
                // unfortunately, the ellipsis sometimes activates too early, meaning that even with clientWidth == scrollWidth some data is hidden by the ellipsis.
                // fortunately, these tooltips are just a convenience so I don't mind if they fail to show.
                // I don't know whether clientWidth or offsetWidth is better. clientWidth was more convenient in Chrome74.
                var e = cell.getElement();
                //return '' + e.offsetWidth + ' || ' + e.scrollWidth + ' || ' + e.clientWidth;
                if (e.clientWidth >= e.scrollWidth) {
                    return false; // all the text is shown, so there is no '...', so no tooltip is needed
                } else {
                    return cell.getValue();
                }
            }
        });
    });
});
