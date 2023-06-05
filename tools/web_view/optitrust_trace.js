// The function [Trace.dump_trace_to_js] in [trace.ml] documents the data representation
// used for the trace. It looks like this, where 0 is always the id of the root step.

// TODO: restore the feature of displaying sources at lines

/*
var startupOpenStep = 45; // optional binding
var steps = [];
steps[0] = {
   id: 0, // corresponds to the step number
   kind: "..",
   exectime: 0.0453;   // in seconds
   name: "..",
   args: [ { name: "..", value: ".."}, { name: "..", value: ".." } ],
   justif: ["..", ".." ],
   isvalid: true,
   script: window.atob("..."),
   scriptline: 23, // possibly undefined
   astBefore: window.atob("..."), // NOT YET IMPLEMENTED; could also an id of an source code stored in a different array, for improved factorization
   astAfter: window.atob("..."), // NOT YET IMPLEMENTED
   diff: window.atob("..."), // could be slow if requested for all!
   sub: [ j1, j2, ... jK ]  // ids of the sub-steps
   }
step[1] = ...

The function initSteps() builds an array of [bigsteps] and an array of [smallsteps],
in order of appearance.

Representation of a big step:
  bigsteps[k] = step object with additional fields:
      id: // the id of the step in the steps array
      start: // an index in smallstep array, inclusive
      stop: // an index in smallstep array, exclusive
  smallsteps[k] = step object with additional fields:
      id: // the id of the step in the steps array

*/

// data structures filled by function initSteps
var smallsteps = [];
var bigsteps = [];
var hasBigsteps = undefined; // false iff bigsteps is empty


// checkbox status; may change default values here
var optionTargetSteps = false;
var optionExectime = false;
var optionIOSteps = false;
var optionShowAST = false;

//---------------------------------------------------
// Code Mirror editor
// Documentation: https://codemirror.net/doc/manual.html

var editor;

// Initialize the editor
function initEditor() {
  editor = CodeMirror.fromTextArea(document.getElementById('source_code'), {
    mode: "text/x-c++src", // or text/x-rustsrc
    lineNumbers: true,
    lineWrapping: true,
    readOnly: true,
    tabSize: 2,
    extraKeys: {
      'N': function(cm) { console.log("pressed N in editor"); },
    },
    // no codemirrror scroll bars, display full code
    scrollbarStyle: null,
    viewportMargin: Infinity,
  });

  // DEPRECATED:
  // Dimensions for codemirror window, if the code has more lines that the
  // size of the window a scrollbar will appear
  // editor.setSize(700, 600);

  // CURRENTLY NOT USED
  // Add a "getSelectedLoc" function
  // Returns a value in the form:
  // { start: { line: 1, col: 4 }, end: { line: 3, col: 8 } }
  editor.getSelectedLoc = function() {
    var from = editor.getCursor(true);
    var to = editor.getCursor(false);
    // Adding 1 because compilers counts from 1, and Codemirror from 0
    return { start: { line: from.line + 1, col: from.ch + 1 },
      end: { line: to.line + 1, col: to.ch + 1 } };
  };
}

function scrollToLoc(loc) { // ({from: from, to: to});
  editor.scrollIntoView(loc, 100);
}

function scrollToLine(line) { // ({from: from, to: to});
  var vMargin = 100; // pixels
  editor.scrollIntoView(line, vMargin);
}

// hook for codemirror
/*
$(document).on('mouseup', '.CodeMirror', function () {
  // triggered on mouseup events in codemirror, eg to update selection
  // editor.getSelectedLoc()
});*/


//---------------------------------------------------
// Diff2html

var configuration = {
   inputFormat: 'json',
   drawFileList: false,
   // fileListToggle: false,
   // fileListStartVisible: false
   fileContentToggle: false,
   matching: 'lines',
   outputFormat: 'side-by-side',
   synchronisedScroll: true,
   highlight: true,
   renderNothingWhenEmpty: false,
   // LATER tune?
   // matchWordsThreshold : similarity threshold for word matching, default is 0.25
   // matchWordsThresholdmatchingMaxComparisons: perform at most this much comparisons for line matching a block of changes, default is 2500
   // maxLineSizeInBlockForComparison: maximum number os characters of the bigger line in a block to apply comparison, default is 200
   // maxLineLengthHighlight: only perform diff changes highlight if lines are smaller than this, default is 10000
   };

// TODO: the horizontal scrollbars should not be needed if the diff contains no long lines.

var hiddenLines = false;
function hideLines() {
  hiddenLines = true;
  $(".d2h-code-side-line").css("padding", "0");
  $(".d2h-code-side-linenumber").css("display", "none");
//   $(".d2h-code-side-linenumber").html("");
//  $(".d2h-code-side-linenumber").css("width", "0em");
}


//---------------------------------------------------

// Could use Underscore's _escape method.
function escapeHTML(s) {
  return s.replaceAll(/&/g,'&amp;')
          .replaceAll(/</g,'&lt;')
          .replaceAll(/>/g,'&gt;');
}

function newlinetobr(s) {
  return s.trim().replaceAll("\n", "<br/>");
}

function htmlButton(id, label, css, onclick) {
  return "<button id='" + id + "' class='" + css + "' type ='button' onclick='" + onclick + "'>" + label + "</button>";
}

function htmlCheckbox(id, label, css, onclick) {
  return "<label><input id='" + id + "' class='" + css + "' type='checkbox' onclick='" + onclick + "'>" + label + "</label>";
}


function htmlSpan(contents, css) {
  return "<span class='" + css +"'>" + contents + "</span>";
}

var maxButtons = 10;

var curSource = -1;
var curSdiff = -1;
var curBdiff = -1;
var idSourceLeft = -1;
var idSourceRight = -1;
var selectedStep = undefined; // stores a step object


function loadDiffFromString(diffString) {
  // this function should be called only after DOM contents is loaded
 var targetElement = document.getElementById("diffDiv");
 var diff2htmlUi = new Diff2HtmlUI(targetElement, diffString, configuration);
 diff2htmlUi.draw();
 diff2htmlUi.highlightCode();

 // identify the two sides of the diff, and register handlers for click on the line numbers;
 $('.d2h-file-side-diff').first().addClass('diffBefore');
 $('.d2h-file-side-diff').last().addClass('diffAfter');
 $('.diffBefore .d2h-code-side-linenumber').click(function(e) {
     var line = e.target.outerText; loadSourceAtLine(selectedStep.ast_before, line); });
 $('.diffAfter .d2h-code-side-linenumber').click(function(e) {
     var line = e.target.outerText; loadSourceAtLine(selectedStep.ast_after, line); });

  // if hideLines() has been called once, call it again
  if (hiddenLines) {
    hideLines();
  }
}

function resetView() {
  $("#sourceDiv").hide();
  $("#diffDiv").hide();
  $("#detailsDiv").html("");
  $("#infoDiv").html("");
  curSource = -1;
  curSdiff = -1;
  curBdiff = -1;
  $(".ctrl-button").removeClass("ctrl-button-selected ctrl-button-covered");
}

function getBdiffCovering(sdiffId) {
  for (var i = 0; i < bigsteps.length; i++) {
    if (bigsteps[i].start <= sdiffId && sdiffId < bigsteps[i].stop) {
      return i;
    }
  }
  return -1; // should not happen
}

function showOrHide(obj, visible) {
  if (visible) {
    obj.show();
  } else {
    obj.hide();
  }
}

function hideNoncoveredButtons(bdiffId) {
  // if (smallsteps.length <= maxButtons)
  //  return; // no hiding needed if only a few steps
  var step = bigsteps[bdiffId];
  var start = step.start;
  var stop = step.stop;
  /*for (var i = 0; i <= codes.length; i++) {
    // DEPRECATED showOrHide($("#button_code_" + i), (start <= i && i <= stop));
    if (i < smallsteps.length)
      showOrHide($("#button_sdiff_" + i), (start <= i && i < stop));
  }*/
  for (var i = 0; i < smallsteps.length; i++) {
    showOrHide($("#button_sdiff_" + i), (start <= i && i < stop));
  }
}

function showBdiffCovered(sdiffId) {
  var bdiffId = getBdiffCovering(sdiffId);
  if (bdiffId == -1)
    return;
  $("#button_bdiff_" + bdiffId).addClass("ctrl-button-covered");
  hideNoncoveredButtons(bdiffId);
}

function displayInfo(descr) {
  $("#infoDiv").html(descr);
}

function loadSourceAtLine(sourceCode, line) {
  loadSource(sourceCode);
  // TODO: autoscrolling does not work...
  var delay = 1000; // milliseconds
   setTimeout(function() { scrollToLine(line); },
   delay);
}

function loadSource(sourceCode, dontResetView) {
  if (! dontResetView) {
    resetView();
  }
  $("#sourceDiv").show();
  editor.setValue(sourceCode);
  // DEPRECTED $("#button_code_" + id).addClass("ctrl-button-selected");
  // showBdiffCovered(id);
  //curSource = id;
}

function loadSdiff(id) {
  $('#button_sdiff_next').focus();
  resetView();
  $("#diffDiv").show();
  var step = smallsteps[id];
  selectedStep = step;
  loadDiffFromString(step.diff);
  var sStep = htmlSpan(newlinetobr(escapeHTML(step.script)), "step-info");
  if (optionExectime) {
    var sTime = htmlSpan(Math.round(1000 * step.exectime) + "ms", "timing-info") + "<div style='clear: both'></div>";
    sStep += sTime;
  }
  displayInfo(sStep);
  $("#button_sdiff_" + id).addClass("ctrl-button-selected");
  // DEPRECATED $("#button_code_" + id).addClass("ctrl-button-covered");
  // DEPRECATED $("#button_code_" + (id+1)).addClass("ctrl-button-covered");
  showBdiffCovered(id);
  curSdiff = id;
  idSourceLeft = id;
  idSourceRight = id+1;
}

function loadBdiff(id) {
  $('#button_bdiff_next').focus();
  resetView();
  $("#diffDiv").show();
  var step = bigsteps[id];
  selectedStep = step;
  loadDiffFromString(step.diff);
  $("#button_bdiff_" + id).addClass("ctrl-button-selected");
  var sStep = htmlSpan(escapeHTML(step.script), "step-info");
  displayInfo(sStep);
  curBdiff = id;
  // $("#button_sdiff_" + bigsteps[id].start).addClass("ctrl-button-covered");
  hideNoncoveredButtons(id);
  curSdiff = bigsteps[id].start - 1; // -1 to anticipate for "next" being pressed
  idSourceLeft = bigsteps[id].start;
  idSourceRight =bigsteps[id].stop;
}

// LATER: simplify, as curSource is deprecated
function nextSdiff() {
  if (curSdiff == -1 && curSource != -1) {
    curSdiff = curSource - 1;
  }
  //var id = Math.min(curSdiff + 1, smallsteps.length-1);
  var id = (curSdiff + 1) % smallsteps.length;
  loadSdiff(id);
}

// LATER: simplify, as curSource is deprecated
function nextBdiff() {
  if (curBdiff == -1) {
    if (curSdiff == -1 && curSource != -1) {
      curSdiff = curSource - 1;
    }
    curBdiff = getBdiffCovering(curSdiff) - 1; // anticipate for the +1 operation
  }
  //var id = Math.min(curBdiff + 1, bigsteps.length-1);
  var id = (curBdiff + 1) % bigsteps.length;
  loadBdiff(id);
}

// handles a click on a step, to view details
function loadStep(idStep) {
  var step = steps[idStep];
  if (optionShowAST) {
    loadSource(step.ast_after, true);
    $("#diffDiv").hide();
    $("#sourceDiv").show();
  } else {
    loadDiffFromString(step.diff);
    $("#diffDiv").show();
    $("#sourceDiv").hide();
  }

  /*
  if (step.kind == "Target") {
    $("#diffDiv").show();
    //loadSource(step.ast_after, true);
  } else {
    $("#diffDiv").hide();
    // loadSource(step.ast_before, true);
    loadSource(step.ast_before, true);
  }
  */
}

function stepToHTML(step) {
  // console.log("steptohtml " + step.id);
  var s = "";
  var sSubs = "";
  for (var i = 0; i < step.sub.length; i++) {
    var substep = steps[step.sub[i]];
    if (!optionTargetSteps && substep.kind == "Target") {
      continue;
    }
    if (!optionIOSteps && substep.kind == "I/O") {
      continue;
    }
    sSubs += "<li>" + stepToHTML(substep) + "</li>\n";
  }
  var validityClass = "";
  if (step.kind == "I/O" || step.kind == "Target") {
    validityClass = "step-io-target";
  } else if (step.kind == "Error") {
    validityClass = "step-error";
  } else if (step.isvalid) {
    validityClass= "step-valid";
  } else {
    validityClass = "step-invalid";
 }
  var sTime = "";
  if (optionExectime) {
    var t = 1000 * step.exectime; // milliseconds
    var nb = "";
    if (t > 10) {
      nb = Math.round(t);
    } else if (t > 1) {
      nb = Math.round(10* t) / 10;
    } else {
      nb = Math.round(100 * t) / 100;
    }
    sTime = "" + nb + "ms";
    var sTimeClass = "exectime-small";
    if (t > 50) {
      sTimeClass = "exectime-heavy";
    } else if (t > 10) {
      sTimeClass = "exectime-mid";
    }
    sTime = "<span class='" + sTimeClass + "'>" + sTime + "</span>";

  }
  var sKind = escapeHTML(step.kind);
  if (step.script_line !== undefined) {
    sKind = "<b>" + step.script_line + "</b>";
  }
  var sScript = escapeHTML(step.script);
  if (step.kind == "Big") {
    sScript = "<b>Bigstep: " + sScript + "</b>";
  }
  var sOnClick = "";
  if (step.hasOwnProperty("id")) { // LATER: refine
    sOnClick = "onclick='loadStep(" + step.id + ")'";
  }

  s += "<div " + sOnClick + " class='step-title " + validityClass + "'>" + sTime + " [" + sKind + "] " + escapeHTML(step.name) + " " + sScript + "</div>";
  for (var i = 0; i < step.justif.length; i++) {
    s += "<div class='step-justif'>" + escapeHTML(step.justif[i]) + "</div>"
  }
  s += "<ul class='step-sub'> " + sSubs + "</ul>\n";
  return s;
}

// handles click on the "details" button
function toggleDetails() {
  var shouldShowDetails = ($("#detailsDiv").html() == "");
  resetView();
  if (shouldShowDetails) {
    $("#diffDiv").hide();
    $("#detailsDiv").html(stepToHTML(selectedStep));
  } else {
    $("#detailsDiv").html("");
  }
}

// handles click on the "all" button
function viewDetailsAll() {
  $("#diffDiv").hide();
  selectedStep = steps[0]; // root
  $("#detailsDiv").html(stepToHTML(selectedStep));
}

function initControls() {
  var s = "";
  function addRow(sTitle, sRow) {
    s += "<span class='row-title'>" + sTitle + ":</span>" + sRow + "<br/>";
  };

  // Code buttons
  /* DEPRECATED
  var sCode = "";
  sCode += htmlButton("button_code_next", "next", "next-button", "nextSource()");
  for (var i = 0; i < codes.length; i++) {
    sCode += htmlButton("button_code_" + i, i, "ctrl-button", "loadSource(" + i + ")");
  }
  addRow("Source", sCode);
  */

  // Big diff buttons
    var sBdiff = "";
  sBdiff += htmlButton("button_bdiff_next", "next", "next-button", "nextBdiff()");
  for (var i = 0; i < bigsteps.length; i++) {
    sBdiff += htmlButton("button_bdiff_" + i, (i+1), "ctrl-button", "loadBdiff(" + i + ")");
  }
  if (bigsteps.length > 1) {
    addRow("BigSteps", sBdiff);
  }

  // Small diff buttons
  var sSdiff = "";
  sSdiff += htmlButton("button_sdiff_next", "next", "next-button", "nextSdiff()");
  for (var i = 0; i < smallsteps.length; i++) {
    sSdiff += htmlButton("button_sdiff_" + i, (i+1), "ctrl-button", "loadSdiff(" + i + ")");
  }
  addRow("SmallSteps", sSdiff);

  // Details button
  s += htmlButton("button_details", "details", "details-button", "toggleDetails()");
  s += htmlButton("button_all", "all", "details-button", "viewDetailsAll()");
  s += htmlCheckbox("option_Exectime", "exectime", "details-checkbox", "updateOptions()");
  s += htmlCheckbox("option_TargetSteps", "target-steps", "details-checkbox", "updateOptions()");
  s += htmlCheckbox("option_IOSteps", "io-steps", "details-checkbox", "updateOptions()");
  s += htmlCheckbox("option_ShowAST", "show-ast", "details-checkbox", "updateOptions()");

  $("#contents").html(s);

  // initialize checkboxes
  $('#option_Exectime').prop('checked', optionExectime);
  $('#option_TargetSteps').prop('checked', optionTargetSteps);
  $('#option_IOSteps').prop('checked', optionIOSteps);
  $('#option_ShowAST').prop('checked', optionShowAST);
}

// handles modification of options by click on the checkboxes
function updateOptions() {
  optionExectime = $('#option_Exectime').prop('checked');
  if (optionExectime) {
    optionIOSteps = true;
    $('#option_IOSteps').prop('checked', optionIOSteps);
  }
  optionTargetSteps = $('#option_TargetSteps').prop('checked');
  optionIOSteps = $('#option_IOSteps').prop('checked');
  optionShowAST = $('#option_ShowAST').prop('checked');
  toggleDetails(); toggleDetails();
}

function initSteps() {
  // reads global variable [steps]
  // writes global variables [smallsteps] and [bigsteps] and [hasBigsteps]
  var rootStep = steps[0];
  var rootSub = rootStep.sub;
  if (rootSub.length == 0) {
    console.log("Error: no steps in tree")
    return;
  }
  // set global variable [hasBigsteps]
  hasBigsteps = (steps[rootSub[0]].kind == "Big");

  // counter used to fill [smallsteps] array
  var curSmallStep = 0;

  // function to iterate over an array of step ids treated as small steps
  function numberSmallSteps(stepIds) {
    for (var i = 0; i < stepIds.length; i++) {
      var smallstep_id = stepIds[i];
      var smallstep = steps[smallstep_id];
      if (smallstep.kind == "I/O") {
        continue;
      }
      if (smallstep.kind != "Small") {
        console.log("Error: numberSmallSteps expected a small step but encountered a step of kind " + smallstep.kind);
      }
      smallstep.id = smallstep_id;
      smallsteps[curSmallStep] = smallstep;
      curSmallStep++;
    }
  };

  if (hasBigsteps) {
    // filling of [bigsteps] and [smallsteps] array
    for (var i = 0; i < rootSub.length; i++) {
      var bigstep_id = rootSub[i];
      var bigstep = steps[bigstep_id];
      if (bigstep.kind == "I/O") {
        continue;
      }
      if (bigstep.kind != "Big") {
        console.log("Error: initSteps expected a big-step but encountered a step of kind " + bigstep.kind);
      }
      bigstep.id = bigstep_id;
      bigstep.start = curSmallStep;
      numberSmallSteps(bigstep.sub);
      bigstep.stop = curSmallStep;
      bigsteps[i] = bigstep;
    }
  } else {
    // there are no big-steps, filling only [smallsteps] array
    numberSmallSteps(rootSub);
  }
}



document.addEventListener('DOMContentLoaded', function () {
  initEditor();
  initSteps();
  initControls();
  initSplitView();
  // editor.setValue("click on a button");
  if (hasBigsteps) {
    loadBdiff(0); }
  else {
    loadSdiff(0);
  }
  $('#button_bdiff_next').focus();

  // start by showing the tree of steps on the root, or the requested step
  var stepInit = 0; // root
  if (typeof startupOpenStep !== "undefined") {
    stepInit = startupOpenStep;
  }
  selectedStep = steps[stepInit];
  toggleDetails(); // calls loadStep(selectedStep)
});

// alternative:
// but there could be many lines..


//// vertical split resizing

var resizer, leftSplit, rightSplit;

function initSplitView() {
  resizer = document.getElementById('vSplitDiv');
  leftSplit = resizer.previousElementSibling;
  rightSplit = resizer.nextElementSibling;
  resizer.addEventListener('mousedown', mouseDownHandler);
}

// The current position of mouse
let x = 0;
let y = 0;

// Width of left side
let leftWidth = 0;

// Handle the mousedown event
// that's triggered when user drags the resizer
const mouseDownHandler = function (e) {
    // Get the current mouse position
    x = e.clientX;
    y = e.clientY;

    leftWidth = leftSplit.getBoundingClientRect().width;

    // Attach the listeners to `document`
    document.addEventListener('mousemove', mouseMoveHandler);
    document.addEventListener('mouseup', mouseUpHandler);
};

const mouseMoveHandler = function (e) {
  // How far the mouse has been moved
  const dx = e.clientX - x;
  const dy = e.clientY - y;

  const newLeftWidth = ((leftWidth + dx) * 100) / resizer.parentNode.getBoundingClientRect().width;
  leftSplit.style.width = `${newLeftWidth}%`;

  document.body.style.cursor = 'col-resize';
  leftSplit.style.userSelect = 'none';
  leftSplit.style.pointerEvents = 'none';

  rightSplit.style.userSelect = 'none';
  rightSplit.style.pointerEvents = 'none';
};

const mouseUpHandler = function () {
  resizer.style.removeProperty('cursor');
  document.body.style.removeProperty('cursor');

  leftSplit.style.removeProperty('user-select');
  leftSplit.style.removeProperty('pointer-events');

  rightSplit.style.removeProperty('user-select');
  rightSplit.style.removeProperty('pointer-events');

  // Remove the handlers of `mousemove` and `mouseup`
  document.removeEventListener('mousemove', mouseMoveHandler);
  document.removeEventListener('mouseup', mouseUpHandler);
};