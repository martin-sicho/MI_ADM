/**
 * Created by sichom on 3/18/16.
 */

initAutocomplete = function(source_url) {
    $(".autocomplete-search").autocomplete({
                source: source_url,
                minLength: 2,
                appendTo: '.autocomplete-search-wrapper'
            }).data("ui-autocomplete")._renderItem = function (ul, item) {

                ul.addClass('collection'); //Ul custom class here

                return $("<li></li>")
                        .addClass('collection-item') //item based custom class to li here
                        .append("<a href='#'>" + item.label + "</a>")
                        .data("ui-autocomplete-item", item)
                        .appendTo(ul);
            };
};