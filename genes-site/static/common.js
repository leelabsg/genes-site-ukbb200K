'use strict';

window._debug = window._debug || {};

// deal with IE11 problems
if (!Math.log10) { Math.log10 = function(x) { return Math.log(x) / Math.LN10; }; }
if (!!window.MSInputMethodContext && !!document.documentMode) { /*ie11*/ $('<style type=text/css>.lz-locuszoom {height: 400px;}</style>').appendTo($('head')); }
if (!String.prototype.includes) {
  String.prototype.includes = function(search, start) {
    'use strict';
    if (typeof start !== 'number') {
      start = 0;
    }
    if (start + search.length > this.length) {
      return false;
    } else {
      return this.indexOf(search, start) !== -1;
    }
  };
}
if (!String.prototype.padStart) {
    String.prototype.padStart = function padStart(targetLength, padString) {
        targetLength = targetLength >> 0; //truncate if number, or convert non-number to 0;
        padString = String(typeof padString !== 'undefined' ? padString : ' ');
        if (this.length >= targetLength) {
            return String(this);
        } else {
            targetLength = targetLength - this.length;
            if (targetLength > padString.length) {
                padString += padString.repeat(targetLength / padString.length); //append to original to ensure we are longer than needed
            }
            return padString.slice(0, targetLength) + String(this);
        }
    };
}
if (!String.prototype.startsWith) {
    Object.defineProperty(String.prototype, 'startsWith', {
        value: function(search, pos) {
            pos = !pos || pos < 0 ? 0 : +pos;
            return this.substring(pos, pos + search.length) === search;
        }
    });
}

(function() {
    // It's unfortunate that these are hard-coded, but it works pretty great, so I won't change it now.
    var autocomplete_bloodhound = new Bloodhound({
        datumTokenizer: Bloodhound.tokenizers.obj.whitespace('display'),
        queryTokenizer: Bloodhound.tokenizers.whitespace,
        identify: function(sugg) { return sugg.display; }, // maybe allows Bloodhound to `.get()`  objects
        remote: {
            url: '/api/autocomplete?query=%QUERY',
            wildcard: '%QUERY',
            rateLimitBy: 'throttle',
            rateLimitWait: 500,
            transform: function(data) {
                // Probably this function reveals that I don't understand Bloodhound.
                // But, I want my previous results to stay around while I keep typing.
                // If the string that's currently in the searchbox matches some string that's been suggested before, I want to see it!
                // This especially happens while I'm typing a chrom-pos-ref-alt.  If what I'm typing agrees with something being suggested, it shouldn't disappear!
                // So, I'm just adding everything to the local index. (Note: NOT localstorage.)
                // Bloodhound appears to perform deduping.
                autocomplete_bloodhound.add(data);
                return data;
            },
        },
        sorter: function(a, b) { return (a.display > b.display) ? 1 : -1; },
    });

    $(function() {
        $('.typeahead').typeahead({
            hint: false,
            highlight: true,
            minLength: 1,
        }, {
            name: 'autocomplete',
            source: autocomplete_bloodhound,
            display: 'value',
            limit: 100,
            templates: {
                suggestion: _.template("<div><%= display %></div>"),
                empty: "<div class='tt-empty-message'>No matches found.</div>"
            }
        });

        $('.typeahead').bind('typeahead:select', function(ev, suggestion) {
            window.location.href = suggestion.url;
        });
    });
})();


// convenience functions
function fmt(format) {
    var args = Array.prototype.slice.call(arguments, 1);
    return format.replace(/{(\d+)}/g, function(match, number) {
        return (typeof args[number] != 'undefined') ? args[number] : match;
    });
}
function deepcopy(obj) {
    return JSON.parse(JSON.stringify(obj));
}
function dataframe_to_objects(df) {
    // convert from dataframe format to javascript's typical array-of-objects format
    // eg, {colA: [10,11], colB: [20,21]} -> [{colA:10, colB:20}, {colA:11, colB: 21}]
    var keys = Object.keys(df);
    var objects = df[keys[0]].map(function() { return {};});
    keys.forEach(function(key) {
        for (var i=0; i<objects.length; i++) {
            objects[i][key] = df[key][i];
        }
    });
    return objects;
}
function objects_to_dataframe(objs) {
    var df = {};
    Object.keys(objs[0]).forEach(function(key) {
        df[key] = objs.map(_.property(key));
    });
    return df;
}

// functions used by many pages
function tabulator_tooltip_maker(cell) {
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
        //return cell.getValue();
        return cell.getElement().innerText;
    }
}

function custom_LocusZoom_Layouts_get(layout_type, layout_name, customizations) {
    // Similar to `LocusZoom.Layouts.get` but also accepts keys like "axes.x.ticks" or "data_layers.1.id_field"
    var layout = LocusZoom.Layouts.get(layout_type, layout_name);
    Object.keys(customizations).forEach(function(key) {
        var value = customizations[key];
        if (!key.includes(".")) {
            layout[key] = value;
        } else {
            var key_parts = key.split(".");
            var obj = layout;
            for (var i=0; i < key_parts.length-1; i++) {
                // TODO: check that `obj` contains `key_parts[i]`
                if (key_parts[i].match(/^[0-9]+$/) && Array.isArray(obj)) {
                    obj = obj[parseInt(key_parts[i])];
                } else {
                    obj = obj[key_parts[i]];
                }
            }
            obj[key_parts[key_parts.length-1]] = value;
        }
    });
    return layout;
}
